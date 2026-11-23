* Run in a fresh Stata session with an older plugin path and this checkout's src directory.
clear all
set more off
version 16.0
args oldplugin package_dir incompatibility
assert `"`oldplugin'"' != "" & `"`package_dir'"' != ""
adopath ++ "`package_dir'"
program _scpc_rust_plugin, plugin using(`"`oldplugin'"')
plugin call _scpc_rust_plugin, ping
if "`incompatibility'" == "methods" {
    plugin call _scpc_rust_plugin, ping_v7
    capture plugin call _scpc_rust_plugin, require_spec 0.3.9
    assert _rc != 0
    local reason "methods spec 0.3.9"
}
else {
    capture plugin call _scpc_rust_plugin, ping_v7
    assert _rc != 0
    local reason "incompatible"
}

input double(s_1 x)
0 -1
1 -4
2  2
4  2
6  1
end
generate double y = x + (s_1 == 0)
generate double z = x + 0.01 * s_1
quietly ivregress 2sls y (x = z), vce(robust)
quietly scpc, avc(.5) k(1) cvs backend(auto)
assert "`e(backend)'" == "fallback"
assert strpos("`e(backend_reason)'", "loaded scpc_rust.plugin") > 0
assert strpos("`e(backend_reason)'", "`reason'") > 0
matrix cv = e(scpccvs)
assert abs(cv[1,3] - 17.432446) < .002
capture scpc, avc(.5) backend(rust)
assert _rc == 198
display "SCPC stale plugin test passed."
