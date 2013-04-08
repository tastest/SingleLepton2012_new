float getUpperLimit_T2tt_LM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 65.0;
     if(seff >= 0.05 && seff < 0.10) ul = 65.7;
     if(seff >= 0.10 && seff < 0.15) ul = 66.8;
     if(seff >= 0.15 && seff < 0.20) ul = 68.2;
     if(seff >= 0.20 && seff < 0.25) ul = 69.5;
     if(seff >= 0.25 && seff < 0.30) ul = 72.1;
     if(seff >= 0.30 && seff < 0.35) ul = 73.5;
     if(seff >= 0.35 && seff < 0.40) ul = 75.6;
     if(seff >= 0.40 && seff < 0.45) ul = 78.5;
     if(seff >= 0.45 && seff < 0.50) ul = 81.5;
     if(seff >= 0.50 && seff < 0.55) ul = 83.7;
     return ul;
}
float getExpectedUpperLimit_T2tt_LM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 75.1;
     if(seff >= 0.05 && seff < 0.10) ul = 75.5;
     if(seff >= 0.10 && seff < 0.15) ul = 76.7;
     if(seff >= 0.15 && seff < 0.20) ul = 78.7;
     if(seff >= 0.20 && seff < 0.25) ul = 80.8;
     if(seff >= 0.25 && seff < 0.30) ul = 83.3;
     if(seff >= 0.30 && seff < 0.35) ul = 86.1;
     if(seff >= 0.35 && seff < 0.40) ul = 89.1;
     if(seff >= 0.40 && seff < 0.45) ul = 93.1;
     if(seff >= 0.45 && seff < 0.50) ul = 96.5;
     if(seff >= 0.50 && seff < 0.55) ul = 99.8;
     return ul;
}
float getExpectedP1UpperLimit_T2tt_LM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 100.6;
     if(seff >= 0.05 && seff < 0.10) ul = 101.3;
     if(seff >= 0.10 && seff < 0.15) ul = 104.5;
     if(seff >= 0.15 && seff < 0.20) ul = 109.7;
     if(seff >= 0.20 && seff < 0.25) ul = 115.1;
     if(seff >= 0.25 && seff < 0.30) ul = 121.3;
     if(seff >= 0.30 && seff < 0.35) ul = 129.0;
     if(seff >= 0.35 && seff < 0.40) ul = 137.0;
     if(seff >= 0.40 && seff < 0.45) ul = 145.0;
     if(seff >= 0.45 && seff < 0.50) ul = 153.7;
     if(seff >= 0.50 && seff < 0.55) ul = 163.4;
     return ul;
}
float getExpectedM1UpperLimit_T2tt_LM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 56.2;
     if(seff >= 0.05 && seff < 0.10) ul = 56.3;
     if(seff >= 0.10 && seff < 0.15) ul = 56.9;
     if(seff >= 0.15 && seff < 0.20) ul = 57.0;
     if(seff >= 0.20 && seff < 0.25) ul = 57.9;
     if(seff >= 0.25 && seff < 0.30) ul = 59.1;
     if(seff >= 0.30 && seff < 0.35) ul = 60.5;
     if(seff >= 0.35 && seff < 0.40) ul = 61.6;
     if(seff >= 0.40 && seff < 0.45) ul = 63.0;
     if(seff >= 0.45 && seff < 0.50) ul = 64.4;
     if(seff >= 0.50 && seff < 0.55) ul = 66.5;
     return ul;
}
