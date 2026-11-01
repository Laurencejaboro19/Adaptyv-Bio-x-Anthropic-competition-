# Block 8: stop a running Modal app (get the id from app list). Used to stop the earlier kras_set4 run.
python -m modal app list
python -m modal app stop ap-XXXXXXXX   # replace with the app id from the list
