typedef struct rvir_input {
  float mass;
  float z;
  float omegam;
  int verbose;
} RVIR_IN;

int parseinput(int argc, char *argv[], RVIR_IN *pinput);
