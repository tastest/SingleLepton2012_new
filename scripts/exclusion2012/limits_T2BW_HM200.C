float getUpperLimit_T2BW_HM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 11.3;
     if(seff >= 0.05 && seff < 0.10) ul = 11.4;
     if(seff >= 0.10 && seff < 0.15) ul = 11.4;
     if(seff >= 0.15 && seff < 0.20) ul = 11.7;
     if(seff >= 0.20 && seff < 0.25) ul = 12.0;
     if(seff >= 0.25 && seff < 0.30) ul = 12.3;
     if(seff >= 0.30 && seff < 0.35) ul = 12.7;
     if(seff >= 0.35 && seff < 0.40) ul = 13.1;
     if(seff >= 0.40 && seff < 0.45) ul = 13.6;
     if(seff >= 0.45 && seff < 0.50) ul = 14.1;
     if(seff >= 0.50 && seff < 0.55) ul = 14.6;
     return ul;
}
float getExpectedUpperLimit_T2BW_HM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 11.8;
     if(seff >= 0.05 && seff < 0.10) ul = 11.9;
     if(seff >= 0.10 && seff < 0.15) ul = 12.0;
     if(seff >= 0.15 && seff < 0.20) ul = 12.3;
     if(seff >= 0.20 && seff < 0.25) ul = 12.5;
     if(seff >= 0.25 && seff < 0.30) ul = 12.9;
     if(seff >= 0.30 && seff < 0.35) ul = 13.3;
     if(seff >= 0.35 && seff < 0.40) ul = 13.8;
     if(seff >= 0.40 && seff < 0.45) ul = 14.2;
     if(seff >= 0.45 && seff < 0.50) ul = 14.8;
     if(seff >= 0.50 && seff < 0.55) ul = 15.4;
     return ul;
}
float getExpectedP1UpperLimit_T2BW_HM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 16.5;
     if(seff >= 0.05 && seff < 0.10) ul = 16.5;
     if(seff >= 0.10 && seff < 0.15) ul = 17.5;
     if(seff >= 0.15 && seff < 0.20) ul = 17.5;
     if(seff >= 0.20 && seff < 0.25) ul = 18.3;
     if(seff >= 0.25 && seff < 0.30) ul = 19.3;
     if(seff >= 0.30 && seff < 0.35) ul = 20.3;
     if(seff >= 0.35 && seff < 0.40) ul = 21.5;
     if(seff >= 0.40 && seff < 0.45) ul = 21.8;
     if(seff >= 0.45 && seff < 0.50) ul = 24.0;
     if(seff >= 0.50 && seff < 0.55) ul = 25.4;
     return ul;
}
float getExpectedM1UpperLimit_T2BW_HM200( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 8.7;
     if(seff >= 0.05 && seff < 0.10) ul = 8.6;
     if(seff >= 0.10 && seff < 0.15) ul = 9.0;
     if(seff >= 0.15 && seff < 0.20) ul = 8.8;
     if(seff >= 0.20 && seff < 0.25) ul = 9.0;
     if(seff >= 0.25 && seff < 0.30) ul = 9.1;
     if(seff >= 0.30 && seff < 0.35) ul = 9.3;
     if(seff >= 0.35 && seff < 0.40) ul = 9.5;
     if(seff >= 0.40 && seff < 0.45) ul = 10.1;
     if(seff >= 0.45 && seff < 0.50) ul = 10.1;
     if(seff >= 0.50 && seff < 0.55) ul = 10.2;
     return ul;
}
