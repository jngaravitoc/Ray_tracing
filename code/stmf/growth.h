#ifndef GROWTH_HINCLUDED
#define GROWTH_HINCLUDED

typedef struct growth_data {
  float omegam;
  float omegal;
  float omegak;
} GROWD;

int growinit(float omegam, float omegal, float omegak, GROWD *pgrowdat);
float dlin0(GROWD growdat);
float dlina(float a, GROWD growdat);
inline float inva3h3(float a, void *growdat);

#endif
