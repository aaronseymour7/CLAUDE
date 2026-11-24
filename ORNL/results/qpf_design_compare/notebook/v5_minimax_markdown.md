### v5: minimax cost (worst-case leakage)

As above, the filter with zero phases acts on the trial state as
 \begin{eqnarray}
 \ket{\Psi_T}\longrightarrow\ket{\phi}&=&\frac{1}{{\cal N}}\prod_i\cos(t_iH)\ket{\Psi_T}=\frac{1}{{\cal N}}\sum_n \bra{n}\Psi_T\rangle\, f(E_n)\ket{n},\qquad f(E)\equiv\prod_i\cos(t_iE)
 \end{eqnarray}
so each eigenstate keeps the fraction $f(E_n)^2$ of its weight, with $f(E_0=0)=1$. For a general trial state with overlap $\gamma=\bra{0}\Psi_T\rangle$ and unknown weights $\lvert\bra{n}\Psi_T\rangle\rvert^2$ on the excited states,
 \begin{eqnarray}
 {\cal N}^2&=&\gamma^2+\sum_{n>0}\lvert\bra{n}\Psi_T\rangle\rvert^2 f(E_n)^2,\\
 1-\lvert\bra{\Psi_{GS}}\phi\rangle\rvert^2&=&\frac{1}{{\cal N}^2}\sum_{n>0}\lvert\bra{n}\Psi_T\rangle\rvert^2 f(E_n)^2
 \;\le\;\frac{1-\gamma^2}{\gamma^2}\,\max_{E\in[\Delta,1]} f(E)^2 .
 \end{eqnarray}
The bound holds for any distribution of the excited weight, and it is reached when that weight sits at the energy where $f^2$ is largest. Instead of assuming a uniform trial state (v4), we minimise the worst leaked level, which needs only the gap $\Delta$:
 \begin{eqnarray}
 C_{v5}&=&\max_{E\in[\Delta,1]}\prod_i\cos^2(t_iE),\qquad \sum_i t_i=T,\quad 0\le t_i\le T/3 .
 \end{eqnarray}
For comparison, v4 minimises
 \begin{eqnarray}
 C_{v4}&=&\left\lvert 1-\frac{1}{\sqrt{\sum_n\prod_i\cos^2(t_iE_n)}}\right\rvert ,
 \end{eqnarray}
which is a monotone function of the mean of $f(E)^2$ over the band (every level weighs the same).

The maximum is not differentiable, so in practice we sample the band on a grid $E_k\in[\Delta,1]$ and minimise its smooth soft-max
 \begin{eqnarray}
 L_\beta&=&\frac{1}{\beta}\log\left[\frac{1}{K}\sum_{k=1}^{K}\exp\!\left(\beta\log f(E_k)^2\right)\right],
 \end{eqnarray}
raising $\beta$ from 2 to 2048. $\beta=1$ gives $\log$ of the mean of $f^2$ (the v4 objective), $\beta\to\infty$ gives $\log\max f^2$ (v5), and $\beta\to0$ gives the mean of $\log f^2$. The final design is checked on a 20 times finer grid.
