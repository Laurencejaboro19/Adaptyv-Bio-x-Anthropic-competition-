# Block 5: BindCraft2's own log for a run, and the files it wrote. Change $tag to look at another run.
$tag = "glp1r_bridge"
cd "$env:USERPROFILE\Downloads\glp1r_modal"
python -m modal run glp1r_app.py::logs --tag $tag
python -m modal run glp1r_app.py::inspect_out --tag $tag
