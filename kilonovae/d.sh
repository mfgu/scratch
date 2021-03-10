#!/bin/bash

python dc.py 60 56
python d.py 60 56 5.0 15.0 10.0 > z60k56.log

python dc.py 60 57
python d.py 60 57 5.0 15.0 10.0 > z60k57.log

python dc.py 60 58
python d.py 60 58 5.0 15.0 10.0 > z60k58.log

python dc.py 60 59
python d.py 60 59 5.0 15.0 5.0 > z60k59.log

python dc.py 60 60
python d.py 60 60 5.0 10.0 2.0 > z60k60.log
