if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YBL088C_1, 5, 20);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_2, 7, 4);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_3, 9, 20);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_4, 7, 15);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_5, 7, 4);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_6, 7, 2);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_7, 9, 20);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_8, 9, 20);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_9, 7, 4);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_10, 9, 15);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_11, 5, 6);
  o += feed_forward(s, nn_yeast_KIN_YBL088C_12, 11, 2);
  o /= 12;
  o = 0.000292749707816895+(0.292722799573086-0.000292749707816895)/(1+exp(100*(0.501975-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYBL088C\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_1, 5, 0);
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_2, 13, 10);
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_3, 5, 6);
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_4, 5, 0);
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_5, 5, 6);
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_6, 5, 2);
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_7, 5, 0);
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_8, 13, 0);
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_9, 5, 0);
  o += feed_forward(s, nn_yeast_KIN_YBR059C_YIL095W_group_10, 5, 6);
  o /= 10;
  o = 5.16581307520249e-07+(0.308144847752964-5.16581307520249e-07)/(1+exp(42.4102*(0.509674-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYBR059C_YIL095W_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YBR136W_1, 7, 10);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_2, 5, 2);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_3, 5, 0);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_4, 7, 2);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_5, 5, 6);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_6, 11, 6);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_7, 11, 2);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_8, 13, 20);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_9, 7, 2);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_10, 5, 2);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_11, 5, 6);
  o += feed_forward(s, nn_yeast_KIN_YBR136W_12, 5, 10);
  o /= 12;
  o = 0.00325819237397371+(0.395458665663462-0.00325819237397371)/(1+exp(73.0682*(0.512668-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYBR136W\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'S') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_1, 7, 20);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_2, 7, 10);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_3, 5, 2);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_4, 7, 10);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_5, 7, 15);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_6, 5, 2);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_7, 11, 20);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_8, 5, 10);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_9, 7, 10);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_10, 7, 10);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_11, 7, 4);
  o += feed_forward(s, nn_yeast_KIN_YHR205W_YBL105C_group_12, 5, 2);
  o /= 12;
  o = 3.99993202370852e-09+(0.129621849014993-3.99993202370852e-09)/(1+exp(100*(0.497422-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYHR205W_YBL105C_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T' || c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_1, 5, 0);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_2, 13, 6);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_3, 11, 15);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_4, 9, 10);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_5, 13, 15);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_6, 11, 15);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_7, 13, 15);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_8, 7, 4);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_9, 11, 10);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_10, 9, 2);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_11, 13, 0);
  o += feed_forward(s, nn_yeast_KIN_YIL035C_YOR061W_group_12, 7, 10);
  o /= 12;
  o = 2.06271259622567e-05+(0.25859204625391-2.06271259622567e-05)/(1+exp(57.443*(0.489669-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYIL035C_YOR061W_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YJL106W_1, 7, 20);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_2, 5, 0);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_3, 5, 15);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_4, 7, 2);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_5, 5, 0);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_6, 7, 4);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_7, 7, 20);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_8, 5, 15);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_9, 7, 15);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_10, 7, 20);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_11, 13, 15);
  o += feed_forward(s, nn_yeast_KIN_YJL106W_12, 5, 2);
  o /= 12;
  o = 1.75545825986181e-07+(0.0687808400529547-1.75545825986181e-07)/(1+exp(99.9998*(0.488864-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYJL106W\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YJR066W_1, 7, 20);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_2, 13, 15);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_3, 11, 6);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_4, 7, 20);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_5, 11, 6);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_6, 7, 10);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_7, 7, 20);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_8, 9, 15);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_9, 7, 10);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_10, 7, 15);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_11, 11, 6);
  o += feed_forward(s, nn_yeast_KIN_YJR066W_12, 7, 10);
  o /= 12;
  o = 0.00678615378657116+(0.0424630896914552-0.00678615378657116)/(1+exp(99.9996*(0.492243-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYJR066W\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YMR001C_1, 11, 2);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_2, 7, 4);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_3, 5, 15);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_4, 13, 4);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_5, 5, 0);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_6, 11, 15);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_7, 13, 4);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_8, 5, 10);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_9, 9, 20);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_10, 11, 15);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_11, 5, 15);
  o += feed_forward(s, nn_yeast_KIN_YMR001C_12, 5, 10);
  o /= 12;
  o = 0.00627737046640504+(0.0538051809128905-0.00627737046640504)/(1+exp(82.0654*(0.484641-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYMR001C\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_1, 9, 6);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_2, 13, 6);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_3, 5, 20);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_4, 11, 10);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_5, 11, 2);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_6, 13, 10);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_7, 9, 6);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_8, 9, 10);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_9, 5, 20);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_10, 7, 10);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_11, 5, 20);
  o += feed_forward(s, nn_yeast_KIN_YPL031C_YBR160W_group_12, 9, 20);
  o /= 12;
  o = 0.0043012645998358+(0.181863344755354-0.0043012645998358)/(1+exp(33.4014*(0.498218-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYPL031C_YBR160W_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_yeast_KIN_YPL153C_1, 5, 15);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_2, 13, 15);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_3, 7, 4);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_4, 5, 4);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_5, 11, 6);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_6, 9, 15);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_7, 13, 6);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_8, 5, 10);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_9, 11, 20);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_10, 5, 10);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_11, 9, 6);
  o += feed_forward(s, nn_yeast_KIN_YPL153C_12, 7, 2);
  o /= 12;
  o = 0.00669217284278663+(0.0914447671110165-0.00669217284278663)/(1+exp(100*(0.490006-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\tyeast\tKIN\tYPL153C\t%.6f\t%.6f\n", o, 0.02);
  }
}

