p = a[0] + x*a[1] + x2*a[2] + x3*a[3] + y*a[4] + x*y*a[5] + x2*y*a[6] + x3*y*a[7] + y2*a[8] + x*y2*a[9] + x2*y2*a[10] + x3*y2*a[11] + y3*a[12] + x*y3*a[13] + x2*y3*a[14] + x3*y3*a[15] + z*a[16] + x*z*a[17] + x2*z*a[18] + x3*z*a[19] + y*z*a[20] + x*y*z*a[21] + x2*y*z*a[22] + x3*y*z*a[23] + y2*z*a[24] + x*y2*z*a[25] + x2*y2*z*a[26] + x3*y2*z*a[27] + y3*z*a[28] + x*y3*z*a[29] + x2*y3*z*a[30] + x3*y3*z*a[31] + z2*a[32] + x*z2*a[33] + x2*z2*a[34] + x3*z2*a[35] + y*z2*a[36] + x*y*z2*a[37] + x2*y*z2*a[38] + x3*y*z2*a[39] + y2*z2*a[40] + x*y2*z2*a[41] + x2*y2*z2*a[42] + x3*y2*z2*a[43] + y3*z2*a[44] + x*y3*z2*a[45] + x2*y3*z2*a[46] + x3*y3*z2*a[47] + z3*a[48] + x*z3*a[49] + x2*z3*a[50] + x3*z3*a[51] + y*z3*a[52] + x*y*z3*a[53] + x2*y*z3*a[54] + x3*y*z3*a[55] + y2*z3*a[56] + x*y2*z3*a[57] + x2*y2*z3*a[58] + x3*y2*z3*a[59] + y3*z3*a[60] + x*y3*z3*a[61] + x2*y3*z3*a[62] + x3*y3*z3*a[63];
px = a[1] + 2*x*a[2] + 3*x2*a[3] + y*a[5] + 2*x*y*a[6] + 3*x2*y*a[7] + y2*a[9] + 2*x*y2*a[10] + 3*x2*y2*a[11] + y3*a[13] + 2*x*y3*a[14] + 3*x2*y3*a[15] + z*a[17] + 2*x*z*a[18] + 3*x2*z*a[19] + y*z*a[21] + 2*x*y*z*a[22] + 3*x2*y*z*a[23] + y2*z*a[25] + 2*x*y2*z*a[26] + 3*x2*y2*z*a[27] + y3*z*a[29] + 2*x*y3*z*a[30] + 3*x2*y3*z*a[31] + z2*a[33] + 2*x*z2*a[34] + 3*x2*z2*a[35] + y*z2*a[37] + 2*x*y*z2*a[38] + 3*x2*y*z2*a[39] + y2*z2*a[41] + 2*x*y2*z2*a[42] + 3*x2*y2*z2*a[43] + y3*z2*a[45] + 2*x*y3*z2*a[46] + 3*x2*y3*z2*a[47] + z3*a[49] + 2*x*z3*a[50] + 3*x2*z3*a[51] + y*z3*a[53] + 2*x*y*z3*a[54] + 3*x2*y*z3*a[55] + y2*z3*a[57] + 2*x*y2*z3*a[58] + 3*x2*y2*z3*a[59] + y3*z3*a[61] + 2*x*y3*z3*a[62] + 3*x2*y3*z3*a[63];
py = a[4] + x*a[5] + x2*a[6] + x3*a[7] + 2*y*a[8] + 2*x*y*a[9] + 2*x2*y*a[10] + 2*x3*y*a[11] + 3*y2*a[12] + 3*x*y2*a[13] + 3*x2*y2*a[14] + 3*x3*y2*a[15] + z*a[20] + x*z*a[21] + x2*z*a[22] + x3*z*a[23] + 2*y*z*a[24] + 2*x*y*z*a[25] + 2*x2*y*z*a[26] + 2*x3*y*z*a[27] + 3*y2*z*a[28] + 3*x*y2*z*a[29] + 3*x2*y2*z*a[30] + 3*x3*y2*z*a[31] + z2*a[36] + x*z2*a[37] + x2*z2*a[38] + x3*z2*a[39] + 2*y*z2*a[40] + 2*x*y*z2*a[41] + 2*x2*y*z2*a[42] + 2*x3*y*z2*a[43] + 3*y2*z2*a[44] + 3*x*y2*z2*a[45] + 3*x2*y2*z2*a[46] + 3*x3*y2*z2*a[47] + z3*a[52] + x*z3*a[53] + x2*z3*a[54] + x3*z3*a[55] + 2*y*z3*a[56] + 2*x*y*z3*a[57] + 2*x2*y*z3*a[58] + 2*x3*y*z3*a[59] + 3*y2*z3*a[60] + 3*x*y2*z3*a[61] + 3*x2*y2*z3*a[62] + 3*x3*y2*z3*a[63];
pz = a[16] + x*a[17] + x2*a[18] + x3*a[19] + y*a[20] + x*y*a[21] + x2*y*a[22] + x3*y*a[23] + y2*a[24] + x*y2*a[25] + x2*y2*a[26] + x3*y2*a[27] + y3*a[28] + x*y3*a[29] + x2*y3*a[30] + x3*y3*a[31] + 2*z*a[32] + 2*x*z*a[33] + 2*x2*z*a[34] + 2*x3*z*a[35] + 2*y*z*a[36] + 2*x*y*z*a[37] + 2*x2*y*z*a[38] + 2*x3*y*z*a[39] + 2*y2*z*a[40] + 2*x*y2*z*a[41] + 2*x2*y2*z*a[42] + 2*x3*y2*z*a[43] + 2*y3*z*a[44] + 2*x*y3*z*a[45] + 2*x2*y3*z*a[46] + 2*x3*y3*z*a[47] + 3*z2*a[48] + 3*x*z2*a[49] + 3*x2*z2*a[50] + 3*x3*z2*a[51] + 3*y*z2*a[52] + 3*x*y*z2*a[53] + 3*x2*y*z2*a[54] + 3*x3*y*z2*a[55] + 3*y2*z2*a[56] + 3*x*y2*z2*a[57] + 3*x2*y2*z2*a[58] + 3*x3*y2*z2*a[59] + 3*y3*z2*a[60] + 3*x*y3*z2*a[61] + 3*x2*y3*z2*a[62] + 3*x3*y3*z2*a[63];