# The 1 comp. PK model with dose-dependent bioavailability.

# [INPUT]
Cl_pop = log(0.2);
V_pop = log(20);
ka_pop = log(0.2);

AGE_med = 69;

beta_V_AGE = 1;
beta_ka_1 = 1;

omega_Cl = 1;
omega_ka = 1;

Cc_a = 1;
Cc_b = 1;

Cl_tv = exp(Cl_pop);
V_tv = exp(V_pop);
ka_tv = exp(ka_pop);


Cl = Cl_tv * exp(omega_Cl)
V = V_tv * (AGE/AGE_med)^beta_V_AGE
ka = ka_tv * exp(omega_ka) * exp(beta_ka_1 * (SEX == 1))

# [MODEL]
# PK model definition

### Priors
MW = 627.57;  					# g/mol for INCB161734

### Explicit functions
Cc = Ac/V;
Cc_nM = Cc/MW*1000000;			# nmol/L

### Initial conditions
Ad(0) = 0;
Ac(0) = 0;

### Differential equations
d/dt(Ad) = - ka*Ad;
d/dt(Ac) = ka*Ad - Cl*Ac/V;
Cc_Res_err = Cc * (1 + Cc_b ) + Cc_a;
