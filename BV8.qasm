OPENQASM 2.0;
include "qelib1.inc";
qreg q[15];
creg c[15];
x q[8];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
cx q[5],q[8];
cx q[6],q[8];
cx q[7],q[8];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];


