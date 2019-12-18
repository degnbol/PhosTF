if (c == 'Y') {
  o = 0;
  o += feed_forward(s, netphosk_InsR_1, 17, 8);
  o += feed_forward(s, netphosk_InsR_2, 9, 8);
  o += feed_forward(s, netphosk_InsR_3, 7, 4);
  o /= 3;
  if (o > 0) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnetphosk\tKIN\tInsR_group\t%.6f\t\n", o);
  }
}
