float getUpperLimit_T2tt_LM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 9.8;
     if(seff >= 0.05 && seff < 0.10) ul = 9.9;
     if(seff >= 0.10 && seff < 0.15) ul = 10.4;
     if(seff >= 0.15 && seff < 0.20) ul = 10.1;
     if(seff >= 0.20 && seff < 0.25) ul = 10.2;
     if(seff >= 0.25 && seff < 0.30) ul = 10.4;
     if(seff >= 0.30 && seff < 0.35) ul = 10.6;
     if(seff >= 0.35 && seff < 0.40) ul = 10.8;
     if(seff >= 0.40 && seff < 0.45) ul = 11.1;
     if(seff >= 0.45 && seff < 0.50) ul = 11.5;
     if(seff >= 0.50 && seff < 0.55) ul = 12.1;
     return ul;
}
float getExpectedUpperLimit_T2tt_LM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 14.5;
     if(seff >= 0.05 && seff < 0.10) ul = 14.6;
     if(seff >= 0.10 && seff < 0.15) ul = 15.2;
     if(seff >= 0.15 && seff < 0.20) ul = 15.1;
     if(seff >= 0.20 && seff < 0.25) ul = 15.4;
     if(seff >= 0.25 && seff < 0.30) ul = 15.8;
     if(seff >= 0.30 && seff < 0.35) ul = 16.2;
     if(seff >= 0.35 && seff < 0.40) ul = 16.8;
     if(seff >= 0.40 && seff < 0.45) ul = 17.4;
     if(seff >= 0.45 && seff < 0.50) ul = 17.9;
     if(seff >= 0.50 && seff < 0.55) ul = 18.9;
     return ul;
}
float getExpectedP1UpperLimit_T2tt_LM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 19.9;
     if(seff >= 0.05 && seff < 0.10) ul = 20.3;
     if(seff >= 0.10 && seff < 0.15) ul = 23.3;
     if(seff >= 0.15 && seff < 0.20) ul = 21.5;
     if(seff >= 0.20 && seff < 0.25) ul = 22.4;
     if(seff >= 0.25 && seff < 0.30) ul = 23.6;
     if(seff >= 0.30 && seff < 0.35) ul = 24.8;
     if(seff >= 0.35 && seff < 0.40) ul = 26.2;
     if(seff >= 0.40 && seff < 0.45) ul = 27.7;
     if(seff >= 0.45 && seff < 0.50) ul = 29.3;
     if(seff >= 0.50 && seff < 0.55) ul = 33.4;
     return ul;
}
float getExpectedM1UpperLimit_T2tt_LM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 10.5;
     if(seff >= 0.05 && seff < 0.10) ul = 10.5;
     if(seff >= 0.10 && seff < 0.15) ul = 11.0;
     if(seff >= 0.15 && seff < 0.20) ul = 10.5;
     if(seff >= 0.20 && seff < 0.25) ul = 10.6;
     if(seff >= 0.25 && seff < 0.30) ul = 10.7;
     if(seff >= 0.30 && seff < 0.35) ul = 10.8;
     if(seff >= 0.35 && seff < 0.40) ul = 10.9;
     if(seff >= 0.40 && seff < 0.45) ul = 11.0;
     if(seff >= 0.45 && seff < 0.50) ul = 11.3;
     if(seff >= 0.50 && seff < 0.55) ul = 11.8;
     return ul;
}
