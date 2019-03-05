qubits 2
map q0, qreg0
map q1, qreg1
map b0, 
map b1, 
PrepZ qreg0
PrepZ qreg1
display
X qreg1
X qreg0
display
H qreg0
display
H qreg1
display
measure qreg0
measure qreg1
display
