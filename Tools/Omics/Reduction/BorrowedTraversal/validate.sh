set -eu
cd /Users/n/numivivo-pca-borrowed-20260912
swiftc -swift-version 6 -O -parse-as-library -I build -I source/Sources/CNumiVivoZlib -I source/Sources/NumiVivoCore/include -L build -lNumiVivoKit Controls.swift build/OmicsHNSW.o build/OmicsGaussian.o build/OmicsMNN.o -framework Accelerate -framework Metal -lc++ -o controls
./controls controls.bin > controls.log
python3 check.py > check.log
