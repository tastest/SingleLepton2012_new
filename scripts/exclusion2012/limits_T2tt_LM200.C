float getUpperLimit_T2tt_LM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 25.3;
     if(seff >= 0.05 && seff < 0.10) ul = 25.7;
     if(seff >= 0.10 && seff < 0.15) ul = 25.7;
     if(seff >= 0.15 && seff < 0.20) ul = 26.2;
     if(seff >= 0.20 && seff < 0.25) ul = 26.6;
     if(seff >= 0.25 && seff < 0.30) ul = 27.5;
     if(seff >= 0.30 && seff < 0.35) ul = 28.2;
     if(seff >= 0.35 && seff < 0.40) ul = 29.2;
     if(seff >= 0.40 && seff < 0.45) ul = 30.1;
     if(seff >= 0.45 && seff < 0.50) ul = 31.0;
     if(seff >= 0.50 && seff < 0.55) ul = 32.0;
     return ul;
}
float getExpectedUpperLimit_T2tt_LM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 31.2;
     if(seff >= 0.05 && seff < 0.10) ul = 31.7;
     if(seff >= 0.10 && seff < 0.15) ul = 32.2;
     if(seff >= 0.15 && seff < 0.20) ul = 32.9;
     if(seff >= 0.20 && seff < 0.25) ul = 33.9;
     if(seff >= 0.25 && seff < 0.30) ul = 34.9;
     if(seff >= 0.30 && seff < 0.35) ul = 36.0;
     if(seff >= 0.35 && seff < 0.40) ul = 37.3;
     if(seff >= 0.40 && seff < 0.45) ul = 38.7;
     if(seff >= 0.45 && seff < 0.50) ul = 40.0;
     if(seff >= 0.50 && seff < 0.55) ul = 41.6;
     return ul;
}
float getExpectedP1UpperLimit_T2tt_LM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 42.2;
     if(seff >= 0.05 && seff < 0.10) ul = 42.7;
     if(seff >= 0.10 && seff < 0.15) ul = 44.0;
     if(seff >= 0.15 && seff < 0.20) ul = 45.6;
     if(seff >= 0.20 && seff < 0.25) ul = 48.0;
     if(seff >= 0.25 && seff < 0.30) ul = 50.5;
     if(seff >= 0.30 && seff < 0.35) ul = 53.5;
     if(seff >= 0.35 && seff < 0.40) ul = 56.8;
     if(seff >= 0.40 && seff < 0.45) ul = 60.2;
     if(seff >= 0.45 && seff < 0.50) ul = 63.9;
     if(seff >= 0.50 && seff < 0.55) ul = 67.4;
     return ul;
}
float getExpectedM1UpperLimit_T2tt_LM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 23.2;
     if(seff >= 0.05 && seff < 0.10) ul = 23.7;
     if(seff >= 0.10 && seff < 0.15) ul = 23.6;
     if(seff >= 0.15 && seff < 0.20) ul = 23.7;
     if(seff >= 0.20 && seff < 0.25) ul = 24.0;
     if(seff >= 0.25 && seff < 0.30) ul = 24.3;
     if(seff >= 0.30 && seff < 0.35) ul = 24.6;
     if(seff >= 0.35 && seff < 0.40) ul = 25.0;
     if(seff >= 0.40 && seff < 0.45) ul = 25.6;
     if(seff >= 0.45 && seff < 0.50) ul = 26.2;
     if(seff >= 0.50 && seff < 0.55) ul = 26.7;
     return ul;
}
