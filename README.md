# CryIoTGS

Data and code repository for Akarsh's FUSars paper on cryogenic TGS during ion irradiation.

## Prerequisites

- Python 3.8 or higher
- Git (for cloning the repository)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/auroraaka/tgs.git
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
   python -m venv tgs-venv
   # On Windows:
   tgs-venv\Scripts\activate
   # On Unix/MacOS:
   source tgs-venv/bin/activate
   
   pip install --upgrade pip
   pip install -e .  # or pip install -r requirements.txt
   ```

## Usage

1. Activate the virtual environment:
   
   **On Unix/MacOS:**
   ```bash
   source tgs-venv/bin/activate
   ```
   
   **On Windows:**
   ```bash
   tgs-venv\Scripts\activate
   ```

2. Edit `config.yaml` with your data path and desired fitting parameters.

3. Run the analysis:
   ```bash
   python main.py
   ```

   Fitting results and figures will be saved in `fit/` and `figures/` directories, respectively.
   You can view example input/output files in the `data/` directory.
