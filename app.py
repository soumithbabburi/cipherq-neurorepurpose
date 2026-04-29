#!/usr/bin/env python3
"""
NeuroRepurpose Intelligence Platform — Entry Point
Run with: python app.py
"""
from dash_app import app, server

if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8501)
