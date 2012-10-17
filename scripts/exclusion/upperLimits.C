float getUpperLimit_SRA( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 183.3;
  if(seff >= 0.05 && seff < 0.15) ul = 187.8;
  if(seff >= 0.15 && seff < 0.25) ul = 195.4;
  if(seff >= 0.25 && seff < 0.35) ul = 207.8;
  if(seff >= 0.35 && seff < 0.45) ul = 222.9;
  if(seff >= 0.45 && seff < 0.55) ul = 234.0;
  return ul;
}


float getExpectedUpperLimit_SRA( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 212.8;
  if(seff >= 0.05 && seff < 0.15) ul = 218.7;
  if(seff >= 0.15 && seff < 0.25) ul = 231.0;
  if(seff >= 0.25 && seff < 0.35) ul = 246.7;
  if(seff >= 0.35 && seff < 0.45) ul = 265.8;
  if(seff >= 0.45 && seff < 0.55) ul = 282.0;
  return ul;
}


float getExpectedP1UpperLimit_SRA( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 284.4;
  if(seff >= 0.05 && seff < 0.15) ul = 295.6;
  if(seff >= 0.15 && seff < 0.25) ul = 323.8;
  if(seff >= 0.25 && seff < 0.35) ul = 365.3;
  if(seff >= 0.35 && seff < 0.45) ul = 411.3;
  if(seff >= 0.45 && seff < 0.55) ul = 515.2;
  return ul;
}


float getExpectedM1UpperLimit_SRA( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 159.3;
  if(seff >= 0.05 && seff < 0.15) ul = 158.6;
  if(seff >= 0.15 && seff < 0.25) ul = 163.9;
  if(seff >= 0.25 && seff < 0.35) ul = 170.3;
  if(seff >= 0.35 && seff < 0.45) ul = 178.9;
  if(seff >= 0.45 && seff < 0.55) ul = 190.5;
  return ul;
}


float getUpperLimit_SRB( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 84.6;
  if(seff >= 0.05 && seff < 0.15) ul = 86.1;
  if(seff >= 0.15 && seff < 0.25) ul = 89.4;
  if(seff >= 0.25 && seff < 0.35) ul = 93.9;
  if(seff >= 0.35 && seff < 0.45) ul = 99.8;
  if(seff >= 0.45 && seff < 0.55) ul = 104.6;
  return ul;
}


float getExpectedUpperLimit_SRB( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 106.8;
  if(seff >= 0.05 && seff < 0.15) ul = 110.2;
  if(seff >= 0.15 && seff < 0.25) ul = 115.2;
  if(seff >= 0.25 && seff < 0.35) ul = 122.4;
  if(seff >= 0.35 && seff < 0.45) ul = 131.2;
  if(seff >= 0.45 && seff < 0.55) ul = 140.7;
  return ul;
}


float getExpectedP1UpperLimit_SRB( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 143.6;
  if(seff >= 0.05 && seff < 0.15) ul = 149.1;
  if(seff >= 0.15 && seff < 0.25) ul = 162.7;
  if(seff >= 0.25 && seff < 0.35) ul = 182.5;
  if(seff >= 0.35 && seff < 0.45) ul = 205.6;
  if(seff >= 0.45 && seff < 0.55) ul = 231.4;
  return ul;
}


float getExpectedM1UpperLimit_SRB( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 79.3;
  if(seff >= 0.05 && seff < 0.15) ul = 78.4;
  if(seff >= 0.15 && seff < 0.25) ul = 80.8;
  if(seff >= 0.25 && seff < 0.35) ul = 83.4;
  if(seff >= 0.35 && seff < 0.45) ul = 87.2;
  if(seff >= 0.45 && seff < 0.55) ul = 90.4;
  return ul;
}

float getUpperLimit_SRC( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 40.7;
  if(seff >= 0.05 && seff < 0.15) ul = 41.7;
  if(seff >= 0.15 && seff < 0.25) ul = 43.9;
  if(seff >= 0.25 && seff < 0.35) ul = 46.6;
  if(seff >= 0.35 && seff < 0.45) ul = 49.8;
  if(seff >= 0.45 && seff < 0.55) ul = 53.3;
  return ul;
}


float getExpectedUpperLimit_SRC( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 46.2;
  if(seff >= 0.05 && seff < 0.15) ul = 47.2;
  if(seff >= 0.15 && seff < 0.25) ul = 49.6;
  if(seff >= 0.25 && seff < 0.35) ul = 53.0;
  if(seff >= 0.35 && seff < 0.45) ul = 56.8;
  if(seff >= 0.45 && seff < 0.55) ul = 61.2;
  return ul;
}


float getExpectedP1UpperLimit_SRC( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 62.5;
  if(seff >= 0.05 && seff < 0.15) ul = 64.8;
  if(seff >= 0.15 && seff < 0.25) ul = 70.8;
  if(seff >= 0.25 && seff < 0.35) ul = 79.2;
  if(seff >= 0.35 && seff < 0.45) ul = 89.1;
  if(seff >= 0.45 && seff < 0.55) ul = 100.3;
  return ul;
}


float getExpectedM1UpperLimit_SRC( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 34.0;
  if(seff >= 0.05 && seff < 0.15) ul = 33.8;
  if(seff >= 0.15 && seff < 0.25) ul = 35.1;
  if(seff >= 0.25 && seff < 0.35) ul = 36.4;
  if(seff >= 0.35 && seff < 0.45) ul = 38.0;
  if(seff >= 0.45 && seff < 0.55) ul = 40.0;
  return ul;
}

float getUpperLimit_SRD( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 28.2;
  if(seff >= 0.05 && seff < 0.15) ul = 29.0;
  if(seff >= 0.15 && seff < 0.25) ul = 30.9;
  if(seff >= 0.25 && seff < 0.35) ul = 33.6;
  if(seff >= 0.35 && seff < 0.45) ul = 36.3;
  if(seff >= 0.45 && seff < 0.55) ul = 39.4;
  return ul;
}


float getExpectedUpperLimit_SRD( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 24.7;
  if(seff >= 0.05 && seff < 0.15) ul = 25.5;
  if(seff >= 0.15 && seff < 0.25) ul = 26.8;
  if(seff >= 0.25 && seff < 0.35) ul = 29.1;
  if(seff >= 0.35 && seff < 0.45) ul = 30.9;
  if(seff >= 0.45 && seff < 0.55) ul = 33.6;
  return ul;
}


float getExpectedP1UpperLimit_SRD( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 34.1;
  if(seff >= 0.05 && seff < 0.15) ul = 35.3;
  if(seff >= 0.15 && seff < 0.25) ul = 38.2;
  if(seff >= 0.25 && seff < 0.35) ul = 44.1;
  if(seff >= 0.35 && seff < 0.45) ul = 47.8;
  if(seff >= 0.45 && seff < 0.55) ul = 54.1;
  return ul;
}


float getExpectedM1UpperLimit_SRD( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 18.3;
  if(seff >= 0.05 && seff < 0.15) ul = 18.6;
  if(seff >= 0.15 && seff < 0.25) ul = 19.5;
  if(seff >= 0.25 && seff < 0.35) ul = 22.4;
  if(seff >= 0.35 && seff < 0.45) ul = 21.4;
  if(seff >= 0.45 && seff < 0.55) ul = 22.6;
  return ul;
}


float getUpperLimit_SRE( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 15.4;
  if(seff >= 0.05 && seff < 0.15) ul = 15.7;
  if(seff >= 0.15 && seff < 0.25) ul = 16.7;
  if(seff >= 0.25 && seff < 0.35) ul = 17.9;
  if(seff >= 0.35 && seff < 0.45) ul = 19.2;
  if(seff >= 0.45 && seff < 0.55) ul = 20.8;
  return ul;
}


float getExpectedUpperLimit_SRE( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 14.5;
  if(seff >= 0.05 && seff < 0.15) ul = 14.7;
  if(seff >= 0.15 && seff < 0.25) ul = 15.6;
  if(seff >= 0.25 && seff < 0.35) ul = 16.7;
  if(seff >= 0.35 && seff < 0.45) ul = 17.7;
  if(seff >= 0.45 && seff < 0.55) ul = 19.0;
  return ul;
}


float getExpectedP1UpperLimit_SRE( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 19.9;
  if(seff >= 0.05 && seff < 0.15) ul = 20.7;
  if(seff >= 0.15 && seff < 0.25) ul = 22.9;
  if(seff >= 0.25 && seff < 0.35) ul = 24.7;
  if(seff >= 0.35 && seff < 0.45) ul = 28.1;
  if(seff >= 0.45 && seff < 0.55) ul = 31.5;
  return ul;
}


float getExpectedM1UpperLimit_SRE( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 10.6;
  if(seff >= 0.05 && seff < 0.15) ul = 11.0;
  if(seff >= 0.15 && seff < 0.25) ul = 12.3;
  if(seff >= 0.25 && seff < 0.35) ul = 12.6;
  if(seff >= 0.35 && seff < 0.45) ul = 12.5;
  if(seff >= 0.45 && seff < 0.55) ul = 12.9;
  return ul;
}

float getUpperLimit_SRF( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 8.2;
  if(seff >= 0.05 && seff < 0.15) ul = 8.3;
  if(seff >= 0.15 && seff < 0.25) ul = 8.6;
  if(seff >= 0.25 && seff < 0.35) ul = 9.1;
  if(seff >= 0.35 && seff < 0.45) ul = 9.8;
  if(seff >= 0.45 && seff < 0.55) ul = 10.4;
  return ul;
}


float getExpectedUpperLimit_SRF( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 8.6;
  if(seff >= 0.05 && seff < 0.15) ul = 8.7;
  if(seff >= 0.15 && seff < 0.25) ul = 9.1;
  if(seff >= 0.25 && seff < 0.35) ul = 9.7;
  if(seff >= 0.35 && seff < 0.45) ul = 10.3;
  if(seff >= 0.45 && seff < 0.55) ul = 11.1;
  return ul;
}


float getExpectedP1UpperLimit_SRF( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 12.1;
  if(seff >= 0.05 && seff < 0.15) ul = 12.3;
  if(seff >= 0.15 && seff < 0.25) ul = 13.5;
  if(seff >= 0.25 && seff < 0.35) ul = 14.8;
  if(seff >= 0.35 && seff < 0.45) ul = 16.1;
  if(seff >= 0.45 && seff < 0.55) ul = 18.7;
  return ul;
}


float getExpectedM1UpperLimit_SRF( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 6.3;
  if(seff >= 0.05 && seff < 0.15) ul = 6.4;
  if(seff >= 0.15 && seff < 0.25) ul = 6.5;
  if(seff >= 0.25 && seff < 0.35) ul = 6.6;
  if(seff >= 0.35 && seff < 0.45) ul = 6.7;
  if(seff >= 0.45 && seff < 0.55) ul = 7.4;
  return ul;
}


float getUpperLimit_SRG( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 4.7;
  if(seff >= 0.05 && seff < 0.15) ul = 4.7;
  if(seff >= 0.15 && seff < 0.25) ul = 4.8;
  if(seff >= 0.25 && seff < 0.35) ul = 5.0;
  if(seff >= 0.35 && seff < 0.45) ul = 5.2;
  if(seff >= 0.45 && seff < 0.55) ul = 5.5;
  return ul;
}


float getExpectedUpperLimit_SRG( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 6.2;
  if(seff >= 0.05 && seff < 0.15) ul = 6.3;
  if(seff >= 0.15 && seff < 0.25) ul = 6.7;
  if(seff >= 0.25 && seff < 0.35) ul = 7.1;
  if(seff >= 0.35 && seff < 0.45) ul = 7.1;
  if(seff >= 0.45 && seff < 0.55) ul = 7.7;
  return ul;
}


float getExpectedP1UpperLimit_SRG( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 8.7;
  if(seff >= 0.05 && seff < 0.15) ul = 8.9;
  if(seff >= 0.15 && seff < 0.25) ul = 10.2;
  if(seff >= 0.25 && seff < 0.35) ul = 10.9;
  if(seff >= 0.35 && seff < 0.45) ul = 11.6;
  if(seff >= 0.45 && seff < 0.55) ul = 13.3;
  return ul;
}


float getExpectedM1UpperLimit_SRG( float seff ){
  float ul = 9999.;

  if(seff < 0.05) ul = 4.6;
  if(seff >= 0.05 && seff < 0.15) ul = 4.6;
  if(seff >= 0.15 && seff < 0.25) ul = 4.6;
  if(seff >= 0.25 && seff < 0.35) ul = 4.7;
  if(seff >= 0.35 && seff < 0.45) ul = 4.8;
  if(seff >= 0.45 && seff < 0.55) ul = 4.9;
  return ul;
}
