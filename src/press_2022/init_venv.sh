python3 -m venv venv
git clone https://github.com/iamgroot42/popcll_torch.git  # Speedup for hamming popcount
mv popcll_torch venv/include
venv/bin/python3 -m pip install --upgrade pip
venv/bin/python3 -m pip install -r requirements.txt --extra-index-url https://download.pytorch.org/whl/cu118
#venv/bin/python3 -m pip install venv/include/popcll_torch/popcll_torch/dist/popcll_torch-1.0-cp39-cp39-linux_x86_64.whl
cd venv/include/popcll_torch/popcll_torch/ && ../../../bin/python3 setup.py install
