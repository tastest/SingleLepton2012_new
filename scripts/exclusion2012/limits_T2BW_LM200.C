float getUpperLimit_T2BW_LM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 34.7;
     if(seff >= 0.05 && seff < 0.10) ul = 35.3;
     if(seff >= 0.10 && seff < 0.15) ul = 35.6;
     if(seff >= 0.15 && seff < 0.20) ul = 35.5;
     if(seff >= 0.20 && seff < 0.25) ul = 36.2;
     if(seff >= 0.25 && seff < 0.30) ul = 37.0;
     if(seff >= 0.30 && seff < 0.35) ul = 37.9;
     if(seff >= 0.35 && seff < 0.40) ul = 38.8;
     if(seff >= 0.40 && seff < 0.45) ul = 40.1;
     if(seff >= 0.45 && seff < 0.50) ul = 41.1;
     if(seff >= 0.50 && seff < 0.55) ul = 42.4;
     return ul;
}
float getExpectedUpperLimit_T2BW_LM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 46.9;
     if(seff >= 0.05 && seff < 0.10) ul = 47.6;
     if(seff >= 0.10 && seff < 0.15) ul = 48.4;
     if(seff >= 0.15 && seff < 0.20) ul = 49.6;
     if(seff >= 0.20 && seff < 0.25) ul = 50.7;
     if(seff >= 0.25 && seff < 0.30) ul = 52.3;
     if(seff >= 0.30 && seff < 0.35) ul = 53.9;
     if(seff >= 0.35 && seff < 0.40) ul = 55.6;
     if(seff >= 0.40 && seff < 0.45) ul = 57.7;
     if(seff >= 0.45 && seff < 0.50) ul = 60.0;
     if(seff >= 0.50 && seff < 0.55) ul = 62.1;
     return ul;
}
float getExpectedP1UpperLimit_T2BW_LM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 64.2;
     if(seff >= 0.05 && seff < 0.10) ul = 64.5;
     if(seff >= 0.10 && seff < 0.15) ul = 66.6;
     if(seff >= 0.15 && seff < 0.20) ul = 69.5;
     if(seff >= 0.20 && seff < 0.25) ul = 72.8;
     if(seff >= 0.25 && seff < 0.30) ul = 76.7;
     if(seff >= 0.30 && seff < 0.35) ul = 81.3;
     if(seff >= 0.35 && seff < 0.40) ul = 86.1;
     if(seff >= 0.40 && seff < 0.45) ul = 91.2;
     if(seff >= 0.45 && seff < 0.50) ul = 96.7;
     if(seff >= 0.50 && seff < 0.55) ul = 102.6;
     return ul;
}
float getExpectedM1UpperLimit_T2BW_LM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 34.8;
     if(seff >= 0.05 && seff < 0.10) ul = 35.4;
     if(seff >= 0.10 && seff < 0.15) ul = 35.6;
     if(seff >= 0.15 && seff < 0.20) ul = 35.4;
     if(seff >= 0.20 && seff < 0.25) ul = 35.5;
     if(seff >= 0.25 && seff < 0.30) ul = 35.8;
     if(seff >= 0.30 && seff < 0.35) ul = 36.5;
     if(seff >= 0.35 && seff < 0.40) ul = 37.0;
     if(seff >= 0.40 && seff < 0.45) ul = 37.9;
     if(seff >= 0.45 && seff < 0.50) ul = 38.5;
     if(seff >= 0.50 && seff < 0.55) ul = 39.2;
     return ul;
}
