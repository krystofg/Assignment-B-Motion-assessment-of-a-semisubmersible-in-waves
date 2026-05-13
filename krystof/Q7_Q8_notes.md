# Q7-Q8 poznamky

Zdroj zadani: `files/Assignment 3 (semi-submersible) F26-2.pdf`

Deadline podle zadani: 13 May 2026, 23:59.

## Relevantni faze

Part A vytvari sea state:

- Q1-Q4: storm maxima and extreme-value statistics.
- Q5-Q6: Bretschneider wave spectrum and spectral moments.

Part B tenhle sea state pouziva pro heave response:

- Q7-Q8 jsou zaklad response modelu: platform mass, added mass a undamped uncoupled heave equation.
- Q9-Q12 potom pouzivaji stejnou effective mass a restoring coefficient pro natural period, RAO, response spectrum a operability statistics.

## Q7 zadani

Otazka: What is the platform mass in tons, and what is the effective total mass including hydrodynamic added mass?

Geometrie z Appendix A:

- T = 18.5 m
- L = 72.3 m
- D = 11.3 m
- a = 30.5 m
- b = 5.40 m
- d = T - b = 13.1 m
- f = T - 2b = 7.70 m

Predpoklady pouzite v MATLAB kodu:

- Pro hydrostaticky vytlak a added mass sloupcu pouzivame delku `f = T - 2b`, tj. sloupce od hladiny po horni hranu pontonu.
- Hloubka `d = T - b` zustava relevantni pro Paullingovu wave-excitation syntézu, protoze je to hloubka k pontoon centreline.
- Pontoons aproximujeme podle double-symmetric top view: dva full-length rectangular boxes s delkou `L` a dva side boxes mezi sloupci s delkou `L - 2D`. Prurez vsech pontonu je `(2b) x (2b)`.
- Hustota morske vody je `rho = 1025 kg/m^3`.

Platform mass:

```text
A_col = pi*(D/2)^2
V_cols = 4*A_col*f
V_pontoons = (2*L + 2*(L - 2D))*(2b)*(2b)
V_disp = V_cols + V_pontoons
M = rho*V_disp
```

Added mass coefficients z Appendix B:

- Columns: circular section, vertical motion, `Cm = 1.00`, `A_R = pi*r_col^2`.
- Pontoons: rectangular/square section, vertical motion, `a/b = 1.00`, `Cm = 1.51`, `A_R = pi*b^2`.

Added mass:

```text
A_33_columns = 4*Cm_col*rho*pi*r_col^2*f
A_33_pontoons = Cm_pont*rho*pi*b^2*(2*L + 2*(L - 2D))
A_33 = A_33_columns + A_33_pontoons
M_eff = M + A_33
```

Ocekavane numericke vysledky z aktualniho MATLAB modelu:

```text
V_disp = 31549.015 m^3
M = 32337.740 tons
A_33 = 37762.247 tons
M_eff = 70099.987 tons
```

## Q8 zadani

Otazka: Write down the equation of uncoupled and undamped heave motion.

Pouzijeme heave displacement `z(t)`, positive upward:

```text
(M + A_33)*z_ddot + C_33*z = F_3(t)
```

kde:

```text
C_33 = rho*g*A_wp
A_wp = 4*pi*(D/2)^2
```

Frequency-domain verze pro harmonic wave excitation:

```text
[C_33 - omega^2*(M + A_33)]*z_hat = F_3_hat(omega)
```

Je to undamped a uncoupled model, takze tam neni `B_33*z_dot` term ani roll/pitch/surge/sway/yaw coupling terms.
