@echo off
cd /d "%~dp0"
python -m pip install -r requirements.txt
python backend_app.py
pause
