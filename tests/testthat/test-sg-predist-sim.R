test_that("sg_predist_sim works ", {
  mod_fin <- RxODE({
    # Doses in mg
    # Time in hours

    ### Parameter values
    # Typical
    Cl_pop = 5;
    V_pop = 180;

    ka_pop = 6;


    # Random effects
    omega_Cl = 0;
    omega_V = 0;
    omega_ka = 0;

    # Residual error
    b = 0;

    ### Parameters
    Cl = Cl_pop * exp(omega_Cl);
    V = V_pop * exp(omega_V);
    ka = ka_pop * exp(omega_ka);

    ### Explicit functions
    Cc = Ac/V;                 # nmol/L

    ### Initial conditions
    Ad(0) = 0;          # mg
    Ac(0) = 0;          # mg

    ### ODEs
    d/dt(Ad) = - ka*Ad;
    d/dt(Ac) = ka*Ad - Cl*Cc ;

    CHECKRUV = b;
    Cc_ResErr = Cc + b*Cc;
  })
  res <- sg_predist_sim(obj1, mod_fin, output = "Cc", npop=100)
  expect_equal(res %>% pull(ID) %>% unique() %>% length(), 100)
  expect_equal(res %>% pull(TIME) %>% unique() %>% length(),
               obj1$SDTAB$TIME %>% unique() %>% length())
})
