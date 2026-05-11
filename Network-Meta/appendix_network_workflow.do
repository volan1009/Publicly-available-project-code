
* Step 1. Configure paths, folders, and logging.
version 15.0
clear all
set more off

local project_dir "v:\PP\Network Meta"
cd "`project_dir'"

capture mkdir figures
capture mkdir results

log using "results\appendix_network_workflow.log", replace text


* Step 2. Check required user-written Stata commands.
capture program drop ensure_network_commands
program define ensure_network_commands
    capture which network
    local has_network = (_rc == 0)

    capture which networkplot
    local has_networkplot = (_rc == 0)

    capture which netfunnel
    local has_netfunnel = (_rc == 0)

    if (`has_network' & `has_networkplot' & `has_netfunnel') {
        global APPENDIX_HAS_NETWORKPLOT 1
        global APPENDIX_HAS_NETFUNNEL 1
        exit
    }

    di as txt "Installing required Stata network meta-analysis commands if needed ..."

    if (!`has_network') {
        capture noisily net install network, from("https://raw.githubusercontent.com/UCL/network/master/package/") replace
        capture which network
        local has_network = (_rc == 0)
    }

    if (!`has_network') {
        capture noisily ssc install network, replace
        capture which network
        local has_network = (_rc == 0)
    }

    if (!`has_network') {
        capture noisily net install st0410_1, from("https://www.stata-journal.com/software/sj15-4") replace
        capture which network
        local has_network = (_rc == 0)
    }

    if (!`has_networkplot' | !`has_netfunnel') {
        capture noisily net install network_graphs, from("https://www.cer-methods.com/Stata") replace
        capture which networkplot
        local has_networkplot = (_rc == 0)
        capture which netfunnel
        local has_netfunnel = (_rc == 0)
    }

    if (!`has_networkplot' | !`has_netfunnel') {
        capture which github
        if (_rc) {
            capture noisily net install github, from("https://haghish.github.io/github/") replace
        }

        capture which github
        if (_rc == 0) {
            capture noisily github install achaimani/network_graphs
            capture program drop networkplot
            capture program drop netfunnel

            capture which networkplot
            local has_networkplot = (_rc == 0)
            capture which netfunnel
            local has_netfunnel = (_rc == 0)
        }
    }

    global APPENDIX_HAS_NETWORKPLOT `has_networkplot'
    global APPENDIX_HAS_NETFUNNEL `has_netfunnel'

    if (!`has_network') {
        di as error "Required command 'network' is still unavailable."
        di as error "Run one of these commands inside Stata, then rerun this do-file:"
        di as error "    net install network, from(https://raw.githubusercontent.com/UCL/network/master/package/) replace"
        di as error "    ssc install network, replace"
        di as error "    net install st0410_1, from(https://www.stata-journal.com/software/sj15-4) replace"
        exit 199
    }

    if (!`has_networkplot') {
        di as txt "networkplot is unavailable; network map export will be skipped."
    }

    if (!`has_netfunnel') {
        di as txt "netfunnel is unavailable; funnel plot export will be skipped."
    }
end


* Step 3. Define a shared schema check.
capture program drop assert_columns
program define assert_columns
    syntax varlist

    foreach var of local varlist {
        capture confirm variable `var'
        if (_rc) {
            di as error "Required column `var' not found."
            exit 111
        }
    }
end


* Step 4. Save graphs in editable and image formats.
capture program drop save_graph_outputs
program define save_graph_outputs
    syntax, BASENAME(string) [WIDTH(integer 2400) HEIGHT(integer 0)]

    graph save "`basename'.gph", replace

    if (`height' > 0) {
        graph export "`basename'.png", replace width(`width') height(`height')
    }
    else {
        graph export "`basename'.png", replace width(`width')
    }
end


* Step 5. Build a project-local netfunnel override with a bottom multi-column legend.
capture program drop ensure_local_netfunnel
program define ensure_local_netfunnel
    syntax, PROJECTDIR(string)

    if ("$APPENDIX_HAS_NETFUNNEL" != "1") {
        exit
    }

    local override_dir "`projectdir'\ado_override"
    local override_file "`override_dir'\netfunnel.ado"

    capture mkdir "`override_dir'"

    capture findfile netfunnel.ado
    if (_rc) {
        di as txt "Skipping the project-local netfunnel override because netfunnel.ado was not found."
        exit
    }

    local source_file "`r(fn)'"

    capture noisily filefilter "`source_file'" "`override_file'", ///
        from("col(6)") ///
        to("position(6) ring(1) col(3) colfirst size(vsmall)") ///
        replace

    if (_rc) {
        di as txt "Using the installed netfunnel because the local legend patch could not be created."
        exit
    }

    if (r(occurrences) == 0) {
        di as txt "Using the installed netfunnel because the expected legend marker was not found in the source file."
        exit
    }

    adopath ++ "`override_dir'"
    capture program drop netfunnel
end


* Step 6. Resize funnel graphs for wider export.
capture program drop format_funnel_graph
program define format_funnel_graph
    capture noisily graph display, xsize(10) ysize(7)

    if (_rc) {
        di as txt "Keeping the default funnel graph size because graph display failed in this Stata setup."
    }
end


* Step 7. Build network and funnel plots for continuous outcomes.
capture program drop run_continuous_network
program define run_continuous_network
    syntax, DATAfile(string) OUTcome(string)

    import delimited "`datafile'", clear
    assert_columns study treatment mean stddev samplesize

    replace treatment = lower(trim(treatment))
    isid study treatment

    network setup mean stddev samplesize, ///
        studyvar(study) ///
        trtvar(treatment) ///
        format(augment) ///
        nocodes

    if ("$APPENDIX_HAS_NETWORKPLOT" == "1") {
        capture noisily network map, improve
        if (_rc == 0) {
            save_graph_outputs, basename("figures\network_`outcome'")
        }
        else {
            di as txt "Skipping network map for `outcome' because network map failed in this Stata setup."
        }
    }
    else {
        di as txt "Skipping network map for `outcome' because networkplot is unavailable."
    }

    if ("$APPENDIX_HAS_NETFUNNEL" == "1") {
        network convert pairs

        netfunnel _y _stderr _t1 _t2, ///
            bycomparison ///
            ytitle("Standard error of mean difference")

        format_funnel_graph
        save_graph_outputs, basename("figures\funnel_`outcome'") width(3200) height(2240)
    }
    else {
        di as txt "Skipping funnel plot for `outcome' because netfunnel is unavailable."
    }
end


* Step 8. Build network and funnel plots for binary outcomes.
capture program drop run_binary_network
program define run_binary_network
    syntax, DATAfile(string) OUTcome(string)

    import delimited "`datafile'", clear
    assert_columns study treatment responders samplesize

    replace treatment = lower(trim(treatment))
    assert responders >= 0
    assert responders <= samplesize
    isid study treatment

    network setup responders samplesize, ///
        studyvar(study) ///
        trtvar(treatment) ///
        format(augment) ///
        nocodes

    if ("$APPENDIX_HAS_NETWORKPLOT" == "1") {
        capture noisily network map, improve
        if (_rc == 0) {
            save_graph_outputs, basename("figures\network_`outcome'")
        }
        else {
            di as txt "Skipping network map for `outcome' because network map failed in this Stata setup."
        }
    }
    else {
        di as txt "Skipping network map for `outcome' because networkplot is unavailable."
    }

    if ("$APPENDIX_HAS_NETFUNNEL" == "1") {
        network convert pairs

        netfunnel _y _stderr _t1 _t2, ///
            bycomparison ///
            ytitle("Standard error of log odds ratio")

        format_funnel_graph
        save_graph_outputs, basename("figures\funnel_`outcome'") width(3200) height(2240)
    }
    else {
        di as txt "Skipping funnel plot for `outcome' because netfunnel is unavailable."
    }
end


* Step 9. Check required input files.
capture confirm file "`project_dir'\ADAS_Cog.csv"
if (_rc) {
    di as error "ADAS_Cog.csv not found in `project_dir'"
    exit 601
}

capture confirm file "`project_dir'\AEs.csv"
if (_rc) {
    di as error "AEs.csv not found in `project_dir'"
    exit 601
}


* Step 10. Run both example analyses.
ensure_network_commands
ensure_local_netfunnel, projectdir("`project_dir'")


run_continuous_network, datafile("`project_dir'\ADAS_Cog.csv") outcome("adas_cog")
run_binary_network, datafile("`project_dir'\AEs.csv") outcome("aes")

log close