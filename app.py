#!/usr/bin/env python3
"""
NeuroRepurpose Intelligence Platform — Entry Point
Run with: python app.py
Opens at: http://localhost:8050
"""
from dash_app import app, server

if __name__ == "__main__":
    app.run(debug=False, host="127.0.0.1", port=8050)
