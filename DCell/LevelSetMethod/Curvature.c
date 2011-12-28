#include "LevelSetMethod.h"

inline double GridFunction2D_CurvLSS( double **p, int i, int j, Coor d )
{
  Coor d2 = {d.x*d.x, d.y*d.y, d.z*d.z};
  double px, py, px2, py2, pxx, pyy, pxy, k;

  pxx = p[-4 + j][-4 + i]/(297.*d2.x) + p[-4 + j][-3 + i]/(1188.*d2.x) - \
(2*p[-4 + j][-2 + i])/(2079.*d2.x) - (17*p[-4 + j][-1 + i])/(8316.*d2.x) - \
(5*p[-4 + j][i])/(2079.*d2.x) - (17*p[-4 + j][1 + i])/(8316.*d2.x) - (2*p[-4 \
+ j][2 + i])/(2079.*d2.x) + p[-4 + j][3 + i]/(1188.*d2.x) + p[-4 + j][4 + \
i]/(297.*d2.x) + p[-3 + j][-4 + i]/(297.*d2.x) + p[-3 + j][-3 + \
i]/(1188.*d2.x) - (2*p[-3 + j][-2 + i])/(2079.*d2.x) - (17*p[-3 + j][-1 + \
i])/(8316.*d2.x) - (5*p[-3 + j][i])/(2079.*d2.x) - (17*p[-3 + j][1 + \
i])/(8316.*d2.x) - (2*p[-3 + j][2 + i])/(2079.*d2.x) + p[-3 + j][3 + \
i]/(1188.*d2.x) + p[-3 + j][4 + i]/(297.*d2.x) + p[-2 + j][-4 + \
i]/(297.*d2.x) + p[-2 + j][-3 + i]/(1188.*d2.x) - (2*p[-2 + j][-2 + \
i])/(2079.*d2.x) - (17*p[-2 + j][-1 + i])/(8316.*d2.x) - (5*p[-2 + \
j][i])/(2079.*d2.x) - (17*p[-2 + j][1 + i])/(8316.*d2.x) - (2*p[-2 + j][2 + \
i])/(2079.*d2.x) + p[-2 + j][3 + i]/(1188.*d2.x) + p[-2 + j][4 + \
i]/(297.*d2.x) + p[-1 + j][-4 + i]/(297.*d2.x) + p[-1 + j][-3 + \
i]/(1188.*d2.x) - (2*p[-1 + j][-2 + i])/(2079.*d2.x) - (17*p[-1 + j][-1 + \
i])/(8316.*d2.x) - (5*p[-1 + j][i])/(2079.*d2.x) - (17*p[-1 + j][1 + \
i])/(8316.*d2.x) - (2*p[-1 + j][2 + i])/(2079.*d2.x) + p[-1 + j][3 + \
i]/(1188.*d2.x) + p[-1 + j][4 + i]/(297.*d2.x) + p[j][-4 + i]/(297.*d2.x) + \
p[j][-3 + i]/(1188.*d2.x) - (2*p[j][-2 + i])/(2079.*d2.x) - (17*p[j][-1 + \
i])/(8316.*d2.x) - (5*p[j][i])/(2079.*d2.x) - (17*p[j][1 + i])/(8316.*d2.x) - \
(2*p[j][2 + i])/(2079.*d2.x) + p[j][3 + i]/(1188.*d2.x) + p[j][4 + \
i]/(297.*d2.x) + p[1 + j][-4 + i]/(297.*d2.x) + p[1 + j][-3 + i]/(1188.*d2.x) \
- (2*p[1 + j][-2 + i])/(2079.*d2.x) - (17*p[1 + j][-1 + i])/(8316.*d2.x) - \
(5*p[1 + j][i])/(2079.*d2.x) - (17*p[1 + j][1 + i])/(8316.*d2.x) - (2*p[1 + \
j][2 + i])/(2079.*d2.x) + p[1 + j][3 + i]/(1188.*d2.x) + p[1 + j][4 + \
i]/(297.*d2.x) + p[2 + j][-4 + i]/(297.*d2.x) + p[2 + j][-3 + i]/(1188.*d2.x) \
- (2*p[2 + j][-2 + i])/(2079.*d2.x) - (17*p[2 + j][-1 + i])/(8316.*d2.x) - \
(5*p[2 + j][i])/(2079.*d2.x) - (17*p[2 + j][1 + i])/(8316.*d2.x) - (2*p[2 + \
j][2 + i])/(2079.*d2.x) + p[2 + j][3 + i]/(1188.*d2.x) + p[2 + j][4 + \
i]/(297.*d2.x) + p[3 + j][-4 + i]/(297.*d2.x) + p[3 + j][-3 + i]/(1188.*d2.x) \
- (2*p[3 + j][-2 + i])/(2079.*d2.x) - (17*p[3 + j][-1 + i])/(8316.*d2.x) - \
(5*p[3 + j][i])/(2079.*d2.x) - (17*p[3 + j][1 + i])/(8316.*d2.x) - (2*p[3 + \
j][2 + i])/(2079.*d2.x) + p[3 + j][3 + i]/(1188.*d2.x) + p[3 + j][4 + \
i]/(297.*d2.x) + p[4 + j][-4 + i]/(297.*d2.x) + p[4 + j][-3 + i]/(1188.*d2.x) \
- (2*p[4 + j][-2 + i])/(2079.*d2.x) - (17*p[4 + j][-1 + i])/(8316.*d2.x) - \
(5*p[4 + j][i])/(2079.*d2.x) - (17*p[4 + j][1 + i])/(8316.*d2.x) - (2*p[4 + \
j][2 + i])/(2079.*d2.x) + p[4 + j][3 + i]/(1188.*d2.x) + p[4 + j][4 + \
i]/(297.*d2.x);

  pyy = p[-4 + j][-4 + i]/(297.*d2.y) + p[-4 + j][-3 + i]/(297.*d2.y) + p[-4 \
+ j][-2 + i]/(297.*d2.y) + p[-4 + j][-1 + i]/(297.*d2.y) + p[-4 + \
j][i]/(297.*d2.y) + p[-4 + j][1 + i]/(297.*d2.y) + p[-4 + j][2 + \
i]/(297.*d2.y) + p[-4 + j][3 + i]/(297.*d2.y) + p[-4 + j][4 + i]/(297.*d2.y) \
+ p[-3 + j][-4 + i]/(1188.*d2.y) + p[-3 + j][-3 + i]/(1188.*d2.y) + p[-3 + \
j][-2 + i]/(1188.*d2.y) + p[-3 + j][-1 + i]/(1188.*d2.y) + p[-3 + \
j][i]/(1188.*d2.y) + p[-3 + j][1 + i]/(1188.*d2.y) + p[-3 + j][2 + \
i]/(1188.*d2.y) + p[-3 + j][3 + i]/(1188.*d2.y) + p[-3 + j][4 + \
i]/(1188.*d2.y) - (2*p[-2 + j][-4 + i])/(2079.*d2.y) - (2*p[-2 + j][-3 + \
i])/(2079.*d2.y) - (2*p[-2 + j][-2 + i])/(2079.*d2.y) - (2*p[-2 + j][-1 + \
i])/(2079.*d2.y) - (2*p[-2 + j][i])/(2079.*d2.y) - (2*p[-2 + j][1 + \
i])/(2079.*d2.y) - (2*p[-2 + j][2 + i])/(2079.*d2.y) - (2*p[-2 + j][3 + \
i])/(2079.*d2.y) - (2*p[-2 + j][4 + i])/(2079.*d2.y) - (17*p[-1 + j][-4 + \
i])/(8316.*d2.y) - (17*p[-1 + j][-3 + i])/(8316.*d2.y) - (17*p[-1 + j][-2 + \
i])/(8316.*d2.y) - (17*p[-1 + j][-1 + i])/(8316.*d2.y) - (17*p[-1 + \
j][i])/(8316.*d2.y) - (17*p[-1 + j][1 + i])/(8316.*d2.y) - (17*p[-1 + j][2 + \
i])/(8316.*d2.y) - (17*p[-1 + j][3 + i])/(8316.*d2.y) - (17*p[-1 + j][4 + \
i])/(8316.*d2.y) - (5*p[j][-4 + i])/(2079.*d2.y) - (5*p[j][-3 + \
i])/(2079.*d2.y) - (5*p[j][-2 + i])/(2079.*d2.y) - (5*p[j][-1 + \
i])/(2079.*d2.y) - (5*p[j][i])/(2079.*d2.y) - (5*p[j][1 + i])/(2079.*d2.y) - \
(5*p[j][2 + i])/(2079.*d2.y) - (5*p[j][3 + i])/(2079.*d2.y) - (5*p[j][4 + \
i])/(2079.*d2.y) - (17*p[1 + j][-4 + i])/(8316.*d2.y) - (17*p[1 + j][-3 + \
i])/(8316.*d2.y) - (17*p[1 + j][-2 + i])/(8316.*d2.y) - (17*p[1 + j][-1 + \
i])/(8316.*d2.y) - (17*p[1 + j][i])/(8316.*d2.y) - (17*p[1 + j][1 + \
i])/(8316.*d2.y) - (17*p[1 + j][2 + i])/(8316.*d2.y) - (17*p[1 + j][3 + \
i])/(8316.*d2.y) - (17*p[1 + j][4 + i])/(8316.*d2.y) - (2*p[2 + j][-4 + \
i])/(2079.*d2.y) - (2*p[2 + j][-3 + i])/(2079.*d2.y) - (2*p[2 + j][-2 + \
i])/(2079.*d2.y) - (2*p[2 + j][-1 + i])/(2079.*d2.y) - (2*p[2 + \
j][i])/(2079.*d2.y) - (2*p[2 + j][1 + i])/(2079.*d2.y) - (2*p[2 + j][2 + \
i])/(2079.*d2.y) - (2*p[2 + j][3 + i])/(2079.*d2.y) - (2*p[2 + j][4 + \
i])/(2079.*d2.y) + p[3 + j][-4 + i]/(1188.*d2.y) + p[3 + j][-3 + \
i]/(1188.*d2.y) + p[3 + j][-2 + i]/(1188.*d2.y) + p[3 + j][-1 + \
i]/(1188.*d2.y) + p[3 + j][i]/(1188.*d2.y) + p[3 + j][1 + i]/(1188.*d2.y) + \
p[3 + j][2 + i]/(1188.*d2.y) + p[3 + j][3 + i]/(1188.*d2.y) + p[3 + j][4 + \
i]/(1188.*d2.y) + p[4 + j][-4 + i]/(297.*d2.y) + p[4 + j][-3 + i]/(297.*d2.y) \
+ p[4 + j][-2 + i]/(297.*d2.y) + p[4 + j][-1 + i]/(297.*d2.y) + p[4 + \
j][i]/(297.*d2.y) + p[4 + j][1 + i]/(297.*d2.y) + p[4 + j][2 + i]/(297.*d2.y) \
+ p[4 + j][3 + i]/(297.*d2.y) + p[4 + j][4 + i]/(297.*d2.y);

  pxy = p[-4 + j][-4 + i]/(225.*d.x*d.y) + p[-4 + j][-3 + i]/(300.*d.x*d.y) + \
p[-4 + j][-2 + i]/(450.*d.x*d.y) + p[-4 + j][-1 + i]/(900.*d.x*d.y) - p[-4 + \
j][1 + i]/(900.*d.x*d.y) - p[-4 + j][2 + i]/(450.*d.x*d.y) - p[-4 + j][3 + \
i]/(300.*d.x*d.y) - p[-4 + j][4 + i]/(225.*d.x*d.y) + p[-3 + j][-4 + \
i]/(300.*d.x*d.y) + p[-3 + j][-3 + i]/(400.*d.x*d.y) + p[-3 + j][-2 + \
i]/(600.*d.x*d.y) + p[-3 + j][-1 + i]/(1200.*d.x*d.y) - p[-3 + j][1 + \
i]/(1200.*d.x*d.y) - p[-3 + j][2 + i]/(600.*d.x*d.y) - p[-3 + j][3 + \
i]/(400.*d.x*d.y) - p[-3 + j][4 + i]/(300.*d.x*d.y) + p[-2 + j][-4 + \
i]/(450.*d.x*d.y) + p[-2 + j][-3 + i]/(600.*d.x*d.y) + p[-2 + j][-2 + \
i]/(900.*d.x*d.y) + p[-2 + j][-1 + i]/(1800.*d.x*d.y) - p[-2 + j][1 + \
i]/(1800.*d.x*d.y) - p[-2 + j][2 + i]/(900.*d.x*d.y) - p[-2 + j][3 + \
i]/(600.*d.x*d.y) - p[-2 + j][4 + i]/(450.*d.x*d.y) + p[-1 + j][-4 + \
i]/(900.*d.x*d.y) + p[-1 + j][-3 + i]/(1200.*d.x*d.y) + p[-1 + j][-2 + \
i]/(1800.*d.x*d.y) + p[-1 + j][-1 + i]/(3600.*d.x*d.y) - p[-1 + j][1 + \
i]/(3600.*d.x*d.y) - p[-1 + j][2 + i]/(1800.*d.x*d.y) - p[-1 + j][3 + \
i]/(1200.*d.x*d.y) - p[-1 + j][4 + i]/(900.*d.x*d.y) - p[1 + j][-4 + \
i]/(900.*d.x*d.y) - p[1 + j][-3 + i]/(1200.*d.x*d.y) - p[1 + j][-2 + \
i]/(1800.*d.x*d.y) - p[1 + j][-1 + i]/(3600.*d.x*d.y) + p[1 + j][1 + \
i]/(3600.*d.x*d.y) + p[1 + j][2 + i]/(1800.*d.x*d.y) + p[1 + j][3 + \
i]/(1200.*d.x*d.y) + p[1 + j][4 + i]/(900.*d.x*d.y) - p[2 + j][-4 + \
i]/(450.*d.x*d.y) - p[2 + j][-3 + i]/(600.*d.x*d.y) - p[2 + j][-2 + \
i]/(900.*d.x*d.y) - p[2 + j][-1 + i]/(1800.*d.x*d.y) + p[2 + j][1 + \
i]/(1800.*d.x*d.y) + p[2 + j][2 + i]/(900.*d.x*d.y) + p[2 + j][3 + \
i]/(600.*d.x*d.y) + p[2 + j][4 + i]/(450.*d.x*d.y) - p[3 + j][-4 + \
i]/(300.*d.x*d.y) - p[3 + j][-3 + i]/(400.*d.x*d.y) - p[3 + j][-2 + \
i]/(600.*d.x*d.y) - p[3 + j][-1 + i]/(1200.*d.x*d.y) + p[3 + j][1 + \
i]/(1200.*d.x*d.y) + p[3 + j][2 + i]/(600.*d.x*d.y) + p[3 + j][3 + \
i]/(400.*d.x*d.y) + p[3 + j][4 + i]/(300.*d.x*d.y) - p[4 + j][-4 + \
i]/(225.*d.x*d.y) - p[4 + j][-3 + i]/(300.*d.x*d.y) - p[4 + j][-2 + \
i]/(450.*d.x*d.y) - p[4 + j][-1 + i]/(900.*d.x*d.y) + p[4 + j][1 + \
i]/(900.*d.x*d.y) + p[4 + j][2 + i]/(450.*d.x*d.y) + p[4 + j][3 + \
i]/(300.*d.x*d.y) + p[4 + j][4 + i]/(225.*d.x*d.y);

  px = -p[-4 + j][-4 + i]/(135.*d.x) - p[-4 + j][-3 + i]/(180.*d.x) - p[-4 + \
j][-2 + i]/(270.*d.x) - p[-4 + j][-1 + i]/(540.*d.x) + p[-4 + j][1 + \
i]/(540.*d.x) + p[-4 + j][2 + i]/(270.*d.x) + p[-4 + j][3 + i]/(180.*d.x) + \
p[-4 + j][4 + i]/(135.*d.x) - p[-3 + j][-4 + i]/(135.*d.x) - p[-3 + j][-3 + \
i]/(180.*d.x) - p[-3 + j][-2 + i]/(270.*d.x) - p[-3 + j][-1 + i]/(540.*d.x) + \
p[-3 + j][1 + i]/(540.*d.x) + p[-3 + j][2 + i]/(270.*d.x) + p[-3 + j][3 + \
i]/(180.*d.x) + p[-3 + j][4 + i]/(135.*d.x) - p[-2 + j][-4 + i]/(135.*d.x) - \
p[-2 + j][-3 + i]/(180.*d.x) - p[-2 + j][-2 + i]/(270.*d.x) - p[-2 + j][-1 + \
i]/(540.*d.x) + p[-2 + j][1 + i]/(540.*d.x) + p[-2 + j][2 + i]/(270.*d.x) + \
p[-2 + j][3 + i]/(180.*d.x) + p[-2 + j][4 + i]/(135.*d.x) - p[-1 + j][-4 + \
i]/(135.*d.x) - p[-1 + j][-3 + i]/(180.*d.x) - p[-1 + j][-2 + i]/(270.*d.x) - \
p[-1 + j][-1 + i]/(540.*d.x) + p[-1 + j][1 + i]/(540.*d.x) + p[-1 + j][2 + \
i]/(270.*d.x) + p[-1 + j][3 + i]/(180.*d.x) + p[-1 + j][4 + i]/(135.*d.x) - \
p[j][-4 + i]/(135.*d.x) - p[j][-3 + i]/(180.*d.x) - p[j][-2 + i]/(270.*d.x) - \
p[j][-1 + i]/(540.*d.x) + p[j][1 + i]/(540.*d.x) + p[j][2 + i]/(270.*d.x) + \
p[j][3 + i]/(180.*d.x) + p[j][4 + i]/(135.*d.x) - p[1 + j][-4 + i]/(135.*d.x) \
- p[1 + j][-3 + i]/(180.*d.x) - p[1 + j][-2 + i]/(270.*d.x) - p[1 + j][-1 + \
i]/(540.*d.x) + p[1 + j][1 + i]/(540.*d.x) + p[1 + j][2 + i]/(270.*d.x) + p[1 \
+ j][3 + i]/(180.*d.x) + p[1 + j][4 + i]/(135.*d.x) - p[2 + j][-4 + \
i]/(135.*d.x) - p[2 + j][-3 + i]/(180.*d.x) - p[2 + j][-2 + i]/(270.*d.x) - \
p[2 + j][-1 + i]/(540.*d.x) + p[2 + j][1 + i]/(540.*d.x) + p[2 + j][2 + \
i]/(270.*d.x) + p[2 + j][3 + i]/(180.*d.x) + p[2 + j][4 + i]/(135.*d.x) - p[3 \
+ j][-4 + i]/(135.*d.x) - p[3 + j][-3 + i]/(180.*d.x) - p[3 + j][-2 + \
i]/(270.*d.x) - p[3 + j][-1 + i]/(540.*d.x) + p[3 + j][1 + i]/(540.*d.x) + \
p[3 + j][2 + i]/(270.*d.x) + p[3 + j][3 + i]/(180.*d.x) + p[3 + j][4 + \
i]/(135.*d.x) - p[4 + j][-4 + i]/(135.*d.x) - p[4 + j][-3 + i]/(180.*d.x) - \
p[4 + j][-2 + i]/(270.*d.x) - p[4 + j][-1 + i]/(540.*d.x) + p[4 + j][1 + \
i]/(540.*d.x) + p[4 + j][2 + i]/(270.*d.x) + p[4 + j][3 + i]/(180.*d.x) + p[4 \
+ j][4 + i]/(135.*d.x);

  py = -p[-4 + j][-4 + i]/(135.*d.y) - p[-4 + j][-3 + i]/(135.*d.y) - p[-4 + \
j][-2 + i]/(135.*d.y) - p[-4 + j][-1 + i]/(135.*d.y) - p[-4 + \
j][i]/(135.*d.y) - p[-4 + j][1 + i]/(135.*d.y) - p[-4 + j][2 + i]/(135.*d.y) \
- p[-4 + j][3 + i]/(135.*d.y) - p[-4 + j][4 + i]/(135.*d.y) - p[-3 + j][-4 + \
i]/(180.*d.y) - p[-3 + j][-3 + i]/(180.*d.y) - p[-3 + j][-2 + i]/(180.*d.y) - \
p[-3 + j][-1 + i]/(180.*d.y) - p[-3 + j][i]/(180.*d.y) - p[-3 + j][1 + \
i]/(180.*d.y) - p[-3 + j][2 + i]/(180.*d.y) - p[-3 + j][3 + i]/(180.*d.y) - \
p[-3 + j][4 + i]/(180.*d.y) - p[-2 + j][-4 + i]/(270.*d.y) - p[-2 + j][-3 + \
i]/(270.*d.y) - p[-2 + j][-2 + i]/(270.*d.y) - p[-2 + j][-1 + i]/(270.*d.y) - \
p[-2 + j][i]/(270.*d.y) - p[-2 + j][1 + i]/(270.*d.y) - p[-2 + j][2 + \
i]/(270.*d.y) - p[-2 + j][3 + i]/(270.*d.y) - p[-2 + j][4 + i]/(270.*d.y) - \
p[-1 + j][-4 + i]/(540.*d.y) - p[-1 + j][-3 + i]/(540.*d.y) - p[-1 + j][-2 + \
i]/(540.*d.y) - p[-1 + j][-1 + i]/(540.*d.y) - p[-1 + j][i]/(540.*d.y) - p[-1 \
+ j][1 + i]/(540.*d.y) - p[-1 + j][2 + i]/(540.*d.y) - p[-1 + j][3 + \
i]/(540.*d.y) - p[-1 + j][4 + i]/(540.*d.y) + p[1 + j][-4 + i]/(540.*d.y) + \
p[1 + j][-3 + i]/(540.*d.y) + p[1 + j][-2 + i]/(540.*d.y) + p[1 + j][-1 + \
i]/(540.*d.y) + p[1 + j][i]/(540.*d.y) + p[1 + j][1 + i]/(540.*d.y) + p[1 + \
j][2 + i]/(540.*d.y) + p[1 + j][3 + i]/(540.*d.y) + p[1 + j][4 + \
i]/(540.*d.y) + p[2 + j][-4 + i]/(270.*d.y) + p[2 + j][-3 + i]/(270.*d.y) + \
p[2 + j][-2 + i]/(270.*d.y) + p[2 + j][-1 + i]/(270.*d.y) + p[2 + \
j][i]/(270.*d.y) + p[2 + j][1 + i]/(270.*d.y) + p[2 + j][2 + i]/(270.*d.y) + \
p[2 + j][3 + i]/(270.*d.y) + p[2 + j][4 + i]/(270.*d.y) + p[3 + j][-4 + \
i]/(180.*d.y) + p[3 + j][-3 + i]/(180.*d.y) + p[3 + j][-2 + i]/(180.*d.y) + \
p[3 + j][-1 + i]/(180.*d.y) + p[3 + j][i]/(180.*d.y) + p[3 + j][1 + \
i]/(180.*d.y) + p[3 + j][2 + i]/(180.*d.y) + p[3 + j][3 + i]/(180.*d.y) + p[3 \
+ j][4 + i]/(180.*d.y) + p[4 + j][-4 + i]/(135.*d.y) + p[4 + j][-3 + \
i]/(135.*d.y) + p[4 + j][-2 + i]/(135.*d.y) + p[4 + j][-1 + i]/(135.*d.y) + \
p[4 + j][i]/(135.*d.y) + p[4 + j][1 + i]/(135.*d.y) + p[4 + j][2 + \
i]/(135.*d.y) + p[4 + j][3 + i]/(135.*d.y) + p[4 + j][4 + i]/(135.*d.y);


  px2 = px * px;
  py2 = py * py;

  k   = (pxx*py2 - 2.*px*py*pxy + pyy*px2) / sqrt((px2+py2)*(px2+py2)*(px2+py2));

  return k;
}
