OPENQASM 2.0;
include "qelib1.inc";
qreg q[15];
creg c[15];
x q[3];
h q[0];
h q[4];
h q[2];
h q[11];
h q[3];
h q[0];
cx q[4],q[3];
cx q[2],q[3];
h q[4];
cx q[11],q[3];
h q[2];
h q[11];
h q[3];
