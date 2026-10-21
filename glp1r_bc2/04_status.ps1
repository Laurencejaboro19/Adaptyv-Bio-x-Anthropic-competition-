# Block 4: progress of the saved runs (any time; works with the PC off while the cloud job runs)
cd "$env:USERPROFILE\Downloads\glp1r_modal"
python -m modal run glp1r_app.py::status
