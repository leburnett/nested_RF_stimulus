"""Entry point for the ephys dashboard."""
from dashboard.app import app
from dashboard.constants import PORT

if __name__ == "__main__":
    app.run(debug=True, port=PORT)
