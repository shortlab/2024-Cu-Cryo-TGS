#!/bin/bash

GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

if [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "win32" ]]; then
    VENV_ACTIVATE="cryiotgs-venv\Scripts\activate"
    PYTHON_CMD="python"
else
    VENV_ACTIVATE="source cryiotgs-venv/bin/activate"
    PYTHON_CMD="python3"
fi

echo -e "${BLUE}Setting up TGS Analysis Tool...${NC}"

if ! command -v $PYTHON_CMD &> /dev/null; then
    echo "Python 3 is required but not installed. Please install Python 3.8 or higher."
    exit 1
fi

echo -e "${BLUE}Creating virtual environment..${NC}"
$PYTHON_CMD -m venv cryiotgs-venv

if [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "win32" ]]; then
    . cryiotgs-venv/Scripts/activate
else
    source cryiotgs-venv/bin/activate
fi

echo -e "${BLUE}Upgrading pip...${NC}"
python -m pip install --upgrade pip

echo -e "${BLUE}Installing dependencies...${NC}"
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
else
    pip install -e .
fi

echo -e "${GREEN}Setup complete! To use the tool:${NC}"
if [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "win32" ]]; then
    echo -e "1. Activate the virtual environment with: ${BLUE}cryiotgs-venv\Scripts\activate${NC}"
else
    echo -e "1. Activate the virtual environment with: ${BLUE}source cryiotgs-venv/bin/activate${NC}"
fi
echo -e "2. Run the codebase with: ${BLUE}python main.py${NC}"
