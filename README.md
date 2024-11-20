# CryIoTGS

Data and code repository for Akarsh's FUSars paper on cryogenic TGS during ion irradiation.

## Prerequisites

- Python 3.8 or higher
- Git (for cloning the repository)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/shortlab/2024-Cu-Cryo-TGS.git
   cd tgs
   ```

2. Run the automated setup script:
   
   **On Unix/MacOS:**
   ```bash
   chmod +x setup.sh
   ./setup.sh
   ```

   **On Windows:**
   ```bash
   bash setup.sh
   ```
   
   Note: On Windows, you'll need Git Bash or a similar bash shell. If you don't have one use the manual installation below.

   **Manual Installation (Windows/Unix/MacOS):**
   ```bash
   python -m venv cryiotgs-venv
   # On Windows:
   cryiotgs-venv\Scripts\activate
   # On Unix/MacOS:
   source cryiotgs-venv/bin/activate
   
   pip install --upgrade pip
   pip install -e .  # or pip install -r requirements.txt
   ```

## Usage

1. Activate the virtual environment:
   
   **On Unix/MacOS:**
   ```bash
   source cryiotgs-venv/bin/activate
   ```
   
   **On Windows:**
   ```bash
   cryiotgs-venv\Scripts\activate
   ```

2. Run the codebase:
   ```bash
   python main.py
   ```
   