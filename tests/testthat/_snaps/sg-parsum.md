# sg-parsum file load works

    Code
      sum_tab
    Output
      # A tibble: 7 x 9
      # Groups:   Parameter [7]
        Parameter Estimate      SE   RSE `95% CI` `CV % (95% CI)` `Shrinkage (var~
        <chr>        <dbl>   <dbl> <dbl> <chr>    <chr>                      <dbl>
      1 ka          0.0715 0.00218  3.05 0.0673 ~ <NA>                       NA   
      2 V          20      0.481    2.4  19.1 - ~ <NA>                       NA   
      3 Cl          0.279  0.0121   4.35 0.255 -~ <NA>                       NA   
      4 omega ka    0.3    0.0219   7.28 0.257 -~ 30.713 (26.2, ~            22.5 
      5 omega V     0.237  0.0176   7.45 0.202 -~ 23.993 (20.4, ~            32.2 
      6 omega Cl    0.434  0.0308   7.1  0.373 -~ 45.477 (38.7, ~             1.73
      7 b           0.0575 0.00112  1.95 0.0553 ~ <NA>                       NA   
      # ... with 2 more variables: OFV <dbl>, AIC <dbl>

---

    Code
      sum_tab
    Output
      # A tibble: 7 x 9
      # Groups:   Parameter [7]
        Parameter Estimate      SE   RSE `95% CI` `CV % (95% CI)` `Shrinkage (var~
        <chr>        <dbl>   <dbl> <dbl> <chr>    <chr>                      <dbl>
      1 ka          0.0715 0.00218  3.05 0.0673 ~ <NA>                       NA   
      2 V          20      0.481    2.4  19.1 - ~ <NA>                       NA   
      3 Cl          0.279  0.0121   4.35 0.255 -~ <NA>                       NA   
      4 omega ka    0.3    0.0219   7.28 0.257 -~ 30.713 (26.2, ~            22.5 
      5 omega V     0.237  0.0176   7.45 0.202 -~ 23.993 (20.4, ~            32.2 
      6 omega Cl    0.434  0.0308   7.1  0.373 -~ 45.477 (38.7, ~             1.73
      7 b           0.0575 0.00112  1.95 0.0553 ~ <NA>                       NA   
      # ... with 2 more variables: OFV <dbl>, AIC <dbl>

