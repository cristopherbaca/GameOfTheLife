## 1. The rules and the program
 
Every cell is alive (1) or dead (0). A generation is computed for all cells at once from the 8 surrounding cells (the **Moore neighbourhood**):
 
| Cell now | Number of live neighbours | Cell next |
|---|---|---|
| alive | 2 or 3 | alive (survival) |
| alive | anything else | dead (loneliness or overcrowding) |
| dead | exactly 3 | alive (birth) |
| dead | anything else | dead |
 
This rule is written **B3/S23** (Birth on 3, Survival on 2 or 3).
 
Controls in your program: `Space` pause/run, `R` random board (paused), `E` erase (paused), left click toggles a cell (paused), and the window can be resized.
 
---
 
## 2. Formal model
 
### 2.1 Board, topology and state space
 
The board has $m$ rows and $n$ columns. Because the code uses `np.roll`, the edges wrap around, so the board is a **discrete torus**
 
$$
\mathbb T_{m,n}=\mathbb Z_m\times\mathbb Z_n ,\qquad (i,j)\sim(i+m,\ j)\sim(i,\ j+n).
$$
 
A **configuration** is a function $x:\mathbb T_{m,n}\to\lbrace 0,1\rbrace$, i.e. an $m\times n$ binary matrix. The set of all configurations is
 
$$
X_{m,n}=\lbrace 0,1\rbrace^{\mathbb T_{m,n}},\qquad |X_{m,n}|=2^{mn}.
$$
 
### 2.2 Local and global maps
 
Let $N_x(c)$ be the number of live neighbours of cell $c=(i,j)$:
 
$$
N_x(i,j)=\sum_{(a,b)\ne(0,0),\ a,b\in\lbrace-1,0,1\rbrace}x(i+a\ \mathrm{mod}\ m,\ j+b\ \mathrm{mod}\ n).
$$
 
The **local rule** $f:\lbrace0,1\rbrace\times\lbrace0,\dots,8\rbrace\to\lbrace0,1\rbrace$ is
 
$$
f(s,N)=\begin{cases}1&\text{if }N=3,\\ 1&\text{if }s=1\text{ and }N=2,\\ 0&\text{otherwise.}\end{cases}
$$
 
The **global map** (one generation) applies $f$ to every cell simultaneously:
 
$$
\Phi:X_{m,n}\to X_{m,n},\qquad \Phi(x)(c)=f\big(x(c), N_x(c)\big).
$$
 
The whole game is the sequence $x,\ \Phi(x),\ \Phi^2(x),\dots$ called the **orbit** (or trajectory) of $x$, where $\Phi^t=\Phi\circ\dots\circ\Phi$ ($t$ times). So the program is a **discrete-time dynamical system** $(X_{m,n},\Phi)$.
 
### 2.3 Dictionary: code and mathematics
 
| Code | Mathematics |
|---|---|
| `board` | configuration $x\in X_{m,n}$ |
| `np.roll(board, d, axis)` | shift operator on the torus (Section 3) |
| `count_neighbors` | $N_x=(K*x)$, a circular convolution with the kernel $K$ minus the centre |
| `next_generation` | global map $\Phi$ |
| `survives \| born` | $\Phi(x)=\mathbb 1[N=3]+x\cdot\mathbb 1[N=2]$ (the two cases are disjoint) |
| `random_board` | sample from a Bernoulli product measure with $\rho=0.25$ |
| `compute_layout` | $n=\lfloor W/6\rfloor,\ m=\lfloor H/6\rfloor$ for a window of $W\times H$ pixels |
| `elapsed_ms`, `interval` | fixed-step time accumulator (Section 12) |
 
---
 
## 3. Counting neighbours: shifts, convolution and Fourier
 
### 3.1 Shifts
 
Let $S_y$ and $S_x$ be the operators that shift a configuration by one row and one column (cyclically). `np.roll(board, 1, axis=0)` is $S_y$, and `np.roll(board, -1, axis=1)` is $S_x^{-1}$. On a torus,
 
$$
S_y^{m}=S_x^{n}=I .
$$
 
The double loop in `count_neighbors` adds $S_y^{a}S_x^{b}$ for the 8 offsets $(a,b)\ne(0,0)$, so
 
$$
\boxed{ N=\Big(\sum_{a,b\in\lbrace-1,0,1\rbrace}S_y^{a}S_x^{b}-I\Big)x }
$$
 
### 3.2 Separability
 
The sum over all 9 offsets factorizes:
 
$$
\sum_{a,b}S_y^{a}S_x^{b}=\big(S_y^{-1}+I+S_y\big)\big(S_x^{-1}+I+S_x\big),
$$
 
so $N=(S_y^{-1}+I+S_y)(S_x^{-1}+I+S_x) x-x$. This means the 9-cell box sum can be computed with **4 shifts** instead of 16 `np.roll` calls (your code rolls twice for each of the 8 neighbours):
 
```python
def count_neighbors_fast(board):
    rows = board + np.roll(board, 1, 0) + np.roll(board, -1, 0)   # 3 cells vertically
    box  = rows  + np.roll(rows,  1, 1) + np.roll(rows,  -1, 1)   # 3x3 box sum
    return box - board                                            # remove the centre
```
 
I checked that it gives exactly the same matrix as `count_neighbors`.
 
### 3.3 Convolution and Fourier diagonalization
 
Equivalently, $N=K*x-x$ with the $3\times3$ all-ones kernel, a **circular convolution**. Circular convolutions are diagonalized by the 2D discrete Fourier transform. The character $e^{2\pi i(ky/m+lx/n)}$ is an eigenvector of both shifts, with eigenvalues $e^{2\pi ik/m}$ and $e^{2\pi il/n}$, hence an eigenvector of the neighbour-count operator with eigenvalue
 
$$
\lambda_{k,l}=\big(1+2\cos\tfrac{2\pi k}{m}\big)\big(1+2\cos\tfrac{2\pi l}{n}\big)-1,\qquad k\in\mathbb Z_m,\ l\in\mathbb Z_n .
$$
 
Consequences:
 
* the largest eigenvalue is $\lambda_{0,0}=3\cdot3-1=8$ (every cell has 8 neighbours),
* the smallest is $\big(3\big)\big(-1\big)-1=-4$, so the **spectrum lies in $[-4, 8]$** (checked numerically: min $-4$, max $8$),
* the eigenvector for $\lambda_{0,0}$ is the constant configuration, which gives a conservation law for the neighbour counts:
$$
\sum_c N_x(c)=8\sum_c x(c),
$$
 
i.e. every live cell is counted as a neighbour exactly 8 times.
 
---
 
## 4. The rule as algebra
 
Because "born" and "survives" can never both be true for the same cell, the rule is a plain sum of two indicator functions:
 
$$
\Phi(x)=\mathbb 1[N_x=3]+x\cdot\mathbb 1[N_x=2].
$$
 
Another equivalent form uses the total over the full $3\times3$ block, $T=N+x$ (the cell itself included):
 
$$
\Phi(x)=1\iff T=3\ \ \text{or}\ \ (T=4\text{ and }x=1).
$$
 
The rule depends on the neighbours only through their **sum**, which makes it **outer-totalistic**; it is also invariant under rotations and reflections of the neighbourhood. Both forms above were checked to be identical to your `survives | born` expression on random boards.
 
---
 
## 5. Determinism: the same matrix always ends the same way
 
**Theorem 1 (determinism).** If two configurations on boards of the same size are equal, $x_0=y_0$, then $\Phi^t(x_0)=\Phi^t(y_0)$ for every $t\ge0$.
 
*Proof.* Induction on $t$. For $t=0$ it is the hypothesis. If $\Phi^t(x_0)=\Phi^t(y_0)$ then $\Phi^{t+1}(x_0)=\Phi\big(\Phi^t(x_0)\big)=\Phi\big(\Phi^t(y_0)\big)=\Phi^{t+1}(y_0)$, because $\Phi$ is a function: each cell's next value depends only on the current values in its $3\times3$ block, with no randomness. $\square$
 
**What this says about your program.**
 
* The *only* random element is the `R` key (`random_board`). The evolution itself is fully deterministic: start twice from the same matrix and you get the same movie, forever. If you call `np.random.seed(42)` before `random_board`, even the random start repeats.
* **The board size is part of the initial condition.** The topology (the torus $\mathbb T_{m,n}$) depends on the window size, and `resize_board` crops or pads the matrix when you resize the window. The same drawing in two windows of different size is **two different dynamical systems**. *Measured*, the R-pentomino (Section 13) placed in the centre of an $s\times s$ torus gives, after 2499 generations:
| Torus size $s\times s$ | Live cells at generation 2499 |
|---|---|
| 40 | 57 |
| 60 | 129 |
| 100 | 145 |
| 200 | 155 |
 
The pattern is identical, the fate is not, because gliders and debris wrap around the edges and collide with other parts of the pattern.
 
* Determinism holds **forward only**. Many different configurations lead to the same next configuration (Section 7), so from a later state you cannot recover which start produced it.
 
---
 
## 6. Finite dynamics: every game ends in a cycle
 
On the torus the state space is finite, $|X_{m,n}|=2^{mn}$. Take any $x_0$ and look at its orbit $x_0,x_1,x_2,\dots$ with $x_{t+1}=\Phi(x_t)$.
 
**Theorem 2 (eventual periodicity).** Every orbit is eventually periodic: there exist a **pre-period** $\tau\ge0$ and a **period** $p\ge1$ with $x_{t+p}=x_t$ for all $t\ge\tau$, and $\tau+p\le2^{mn}$.
 
*Proof.* Among $x_0,\dots,x_{2^{mn}}$ there are $2^{mn}+1$ configurations in a set of size $2^{mn}$, so two are equal by the pigeonhole principle: $x_a=x_b$ with $a<b$. By Theorem 1 the orbit repeats from $a$ on with period $b-a$. $\square$
 
In graph terms, $\Phi$ defines a **functional graph** on $X_{m,n}$ in which every node has exactly one outgoing edge. Every connected component consists of a single cycle (an **attractor**) with trees of transient states feeding into it. The cycle lengths are the periods of the "final fates": period 1 is a still life (a stable board, possibly the empty one), period $p>1$ is a blinking pattern, and larger periods with moving parts are the rotating or wrapping configurations of Section 9.
 
### 6.1 How large is the state space?
 
The number of decimal digits of $2^{mn}$ is $\lfloor mn\log_{10}2\rfloor+1$. For a $1920\times1080$ window with `CELL_SIZE = 6` the board is $180\times320=57{,}600$ cells, so
 
$$
|X|=2^{57600}\approx10^{17339}.
$$
 
The theoretical bound $\tau+p\le2^{57600}$ is useless in practice, so the real behaviour has to be explored experimentally.
 
### 6.2 Exact computation on tiny tori
 
For small boards the whole functional graph can be built by brute force (code in Section 14). *Measured*:
 
| Torus | States $2^{mn}$ | States with no predecessor | Cyclic states | Number of attractors | Longest transient |
|---|---|---|---|---|---|
| $4\times4$ | 65,536 | 47,657 (72.7%) | 1,245 | 345 | 9 steps |
| $4\times5$ | 1,048,576 | 769,411 (73.4%) | 1,097 | 461 | 21 steps |
| $5\times5$ | 33,554,432 | 25,033,436 (74.6%) | 8,406 | 5,181 | 51 steps |
 
Cycle structure (period: how many distinct cycles have that period):
 
* $4\times4$: period 1: 53, period 2: 180, period 4: 16, period 8: 96.
* $4\times5$: period 1: 313, 2: 32, 4: 20, 5: 8, 6: 60, 8: 20, 10: 8.
* $5\times5$: period 1: 3456, 2: 1225, 3: 200, 4: 200, 5: 20, 10: 60, 20: 20.
Three things stand out: almost all states are *transient* (only a tiny fraction lies on cycles), roughly three quarters of all states cannot be reached from anything (Section 7), and the maximum transient length grows fast with the board size.
 
---
 
## 7. Irreversibility and Gardens of Eden
 
### 7.1 $\Phi$ is not injective
 
Two different configurations can have the same successor. Example: the empty board and a board with a single live cell both go to the empty board in one step, since an isolated cell has 0 neighbours and dies:
 
$$
\Phi(\text{single cell})=\Phi(\text{empty})=\text{empty}.
$$
 
So $\Phi$ is **not invertible**: the game has an arrow of time. Information is destroyed, and the entropy of the configuration set can only shrink: $|\Phi(X)|<|X|$. For the tori above the image is exactly the complement of the "no predecessor" column (e.g. $17{,}879$ of $65{,}536$ states for $4\times4$).
 
### 7.2 Gardens of Eden
 
A configuration with **no predecessor**, i.e. $x\notin\Phi(X)$, is called a **Garden of Eden**: it can exist only as a starting position, never as the result of a generation.
 
**Corollary (finite board).** On any torus $\mathbb T_{m,n}$, since $\Phi$ is not injective and $X$ is finite, $\Phi$ is not surjective, so Gardens of Eden exist. Their fraction was measured above (about 73-75% on small tori).
 
**Moore–Myhill theorem (infinite plane).** For a cellular automaton on $\mathbb Z^2$, the global map is surjective if and only if it is **pre-injective**, meaning that any two configurations that differ in only finitely many cells and have the same image must be equal. Our example above (empty vs single cell) differs in one cell and has the same image, so $\Phi$ is **not pre-injective**; therefore it is **not surjective**, and Gardens of Eden exist on the infinite plane as well.
 
Every attempt to "run the game backwards" is therefore ambiguous (many predecessors) or impossible (none), and deciding even **surjectivity** is undecidable for general two-dimensional cellular automata (Kari, 1994).
 
---
 
## 8. Change something and everything changes
 
This section covers your question: if you change one cell, what happens later?
 
### 8.1 Locality: nothing travels faster than one cell per generation
 
Define the **Chebyshev distance** $d_\infty(c,c')=\max(|i-i'|,|j-j'|)$ (measured with wrap-around on the torus), and the "ball" $B_t(D)=\lbrace c:\ d_\infty(c,D)\le t\rbrace$.
 
**Theorem 3 (light cone).** If configurations $x$ and $y$ differ only on a set $D$ of cells, then $\Phi^t(x)$ and $\Phi^t(y)$ differ only on $B_t(D)$.
 
*Proof.* For $t=1$: the value $\Phi(x)(c)$ depends only on the $3\times3$ block around $c$. If that block does not meet $D$, then $x$ and $y$ agree on it and so do $\Phi(x)(c)$ and $\Phi(y)(c)$. Hence the difference lies in $B_1(D)$. Induct on $t$, applying the case $t=1$ to $\Phi^{t}(x),\Phi^{t}(y)$. $\square$
 
So there is a maximal "speed of light" $c=1$ cell per generation. For a single flipped cell, the number of differing cells after $t$ generations is at most $(2t+1)^2$ (until the torus is covered).
 
I verified this empirically: in 50 random soups with one cell flipped, the set of differing cells always stayed inside the square of radius $t$ around the flipped cell, for $t=1,\dots,40$.
 
### 8.2 A topological statement: Curtis–Hedlund–Lyndon
 
On the infinite plane, put the product topology on $\lbrace0,1\rbrace^{\mathbb Z^2}$, equivalently the metric
 
$$
d(x,y)=2^{-\min\lbrace \lVert c\rVert_\infty :\ x(c)\ne y(c)\rbrace}.
$$
 
Two configurations are close when they agree on a large square around the origin. The **Curtis–Hedlund–Lyndon theorem** says that a map is the global map of a cellular automaton **if and only if** it is continuous and commutes with all translations. This is the formal statement of Theorem 3: close inputs have close outputs *for any finite number of steps*.
 
But "continuous" does not mean "stable". For the infinite time horizon a tiny change can still grow without bound, because continuity only controls a **bounded region over a bounded time**.
 
### 8.3 The rule is not linear, so there is no superposition
 
Some cellular automata are **linear over $\mathbb F_2$** (for example the "XOR of neighbours" rule). For them, $\Phi(x\oplus y)=\Phi(x)\oplus\Phi(y)$, so the difference between two trajectories is itself a trajectory and the effect of a change can be computed separately from everything else.
 
The Game of Life is **not** linear. Counter-example (found by search; rows shown inside a $5\times5$ window): let $x$ be two vertical live cells, and $y$ a single live cell away from them.
 
```
x          y          x xor y        Φ(x xor y)     Φ(x) xor Φ(y)
.....      .....      .....          .....          .....
.....      .O...      .O...          .....          .....
...O.      .....      ...O.          ..O..          .....
...O.      .....      ...O.          .....          .....
.....      .....      .....          .....          .....
```
 
Each of $x$ and $y$ dies completely in one step, so $\Phi(x)\oplus\Phi(y)$ is empty. But their union has a cell with exactly 3 neighbours, which is **born**. Therefore $\Phi(x\oplus y)\ne\Phi(x)\oplus\Phi(y)$: the effect of a change depends on the context, and a change cannot be analysed in isolation. This is the mathematical root of the "butterfly effect" here.
 
### 8.4 Measured sensitivity on random soups
 
Experiment: 300 random soups with density $0.25$ on a $100\times100$ torus; flip exactly one random cell in each; compare with the unperturbed run by the **Hamming distance** $H_t=\#\lbrace c:\ \Phi^t(x)(c)\ne\Phi^t(y)(c)\rbrace$.
 
| Generation $t$ | Runs where the two coincide again | Median $H_t$ | Mean $H_t$ |
|---|---|---|---|
| 0 | 0% | 1 | 1.0 |
| 1 | 2.3% | 4 | 4.3 |
| 5 | 3.3% | 25 | 24.0 |
| 10 | 3.3% | 63 | 62.1 |
| 25 | 3.3% | 208.5 | 204.2 |
| 50 | 3.3% | 453 | 448.9 |
| 100 | 3.3% | 971 | 938.0 |
| 200 | 3.3% | 1271 | 1232.4 |
| 400 | 3.3% | 1079.5 | 1050.4 |
| 600 | 3.3% | 953 | 930.3 |
 
What the numbers show:
 
* **Healing is rare but real.** In 10 of 300 runs (3.3%) the change disappeared within 3 generations (a typical mechanism: a block missing one cell becomes an L-shape of 3 cells, and the missing corner has exactly 3 neighbours, so it is born and the block is restored in one generation). Once the two orbits coincide they stay equal forever, by Theorem 1.
* **Otherwise the damage spreads** roughly as a growing region, always within the bound $(2t+1)^2$ (for example at $t=10$ the bound is 441 and the median is 63).
* **It then saturates** near a value consistent with *two independent soups*: with final density $q\approx0.05$, independent configurations differ in about $2q(1-q) mn\approx950$ cells, matching the measured $\approx930$-$950$. After a few hundred generations the single-cell change has effectively produced an unrelated evolution.
So the intuition "change one thing and everything changes later" is correct for typical dense starts, with the precise mathematical content: the influence region grows at speed $\le1$, nothing is additive, and in the long run the two histories become statistically unrelated, while the rare exceptions (healing) are permanent.
 
---
 
## 9. Symmetries, spaceships and the torus
 
### 9.1 Symmetries commute with the dynamics
 
Let $g$ be a translation of the torus, or a rotation by $90^\circ$, or a reflection (the symmetry group of the square lattice, $D_4$, together with translations). Because $f$ treats the 8 neighbours identically, 
 
$$
\Phi(g\cdot x)=g\cdot\Phi(x).
$$
 
**Corollary (symmetry is preserved).** If $g\cdot x=x$ then $g\cdot\Phi^t(x)=\Phi^t(g\cdot x)=\Phi^t(x)$. A starting pattern with a symmetry keeps it forever, *no matter how chaotic the evolution looks*. I checked it: a random board made symmetric under all 8 symmetries of the square stayed exactly symmetric for 60 generations. This is why patterns such as the pulsar look so regular. It also means that to break symmetry you must change the starting pattern asymmetrically.
 
### 9.2 Oscillators and spaceships as fixed points of $\Phi^p$
 
Let $\tau_v$ denote translation by a vector $v$.
 
* A **still life** satisfies $\Phi(x)=x$.
* An **oscillator of period $p$** satisfies $\Phi^p(x)=x$ (and no smaller $p$ works).
* A **spaceship** satisfies $\Phi^p(x)=\tau_v x$ with $v\ne0$: it reappears shifted by $v$ after $p$ generations. Its speed is $|v|/p$ cells per generation (in the Chebyshev norm, compared with the light speed $c=1$).
Examples (*measured*): the **glider** has $p=4$, $v=(1,1)$, so speed $c/4$ along the diagonal; the **lightweight spaceship** has $p=4$, $|v|=2$ along a row, speed $c/2$. It is a known result that no orthogonal spaceship can exceed $c/2$.
 
### 9.3 Spaceships on a torus come back
 
In your program the glider does not escape: it wraps around. On $\mathbb T_{m,n}$ the shift by $(1,1)$ per 4 generations returns to the starting position after $t$ generations exactly when $t/4$ is a multiple of both $m$ and $n$:
 
$$
\boxed{ T_{\text{glider}}=4\cdot\mathrm{lcm}(m,n) }
$$
 
Verified: $20\times30\to240$, $15\times15\to60$, $12\times18\to144$, $7\times9\to252$ (all equal to $4 \mathrm{lcm}(m,n)$). If $m$ and $n$ are coprime the glider needs $4mn$ generations, the longest possible wait. A lone glider on a torus is therefore an orbit of exact period $4 \mathrm{lcm}(m,n)$; with a second object in the way, the interaction usually creates something new.
 
---
 
## 10. Random starts: density dynamics
 
### 10.1 What `random_board` samples
 
Each of the $mn$ cells is independently alive with probability $\rho=0.25$ (a Bernoulli product measure). The initial population $P_0\sim\mathrm{Binomial}(mn,\rho)$ satisfies
 
$$
\mathbb E[P_0]=\rho mn,\qquad \mathrm{Var}[P_0]=mn \rho(1-\rho).
$$
 
For the $180\times320$ board: mean $14{,}400$ cells, standard deviation $\sqrt{57600\cdot0.25\cdot0.75}\approx104$ cells, a relative fluctuation of about 0.7%. All random boards therefore start with nearly the same density.
 
### 10.2 The first generation is exact
 
For one generation the neighbours of a cell are 8 independent Bernoulli($p$) variables, so $N\sim\mathrm{Binomial}(8,p)$ and
 
$$
p'=\Pr[\text{alive next}]=\underbrace{\Pr[N=3]}_{\binom83p^3(1-p)^5}+p\cdot\underbrace{\Pr[N=2]}_{\binom82p^2(1-p)^6}.
$$
 
Writing out the binomial coefficients $\binom83=56$ and $\binom82=28$ and combining terms:
 
$$
\boxed{ F(p)=84 p^{3}(1-p)^{6}+56 p^{4}(1-p)^{5} }
$$
 
so $F(0.25)=0.2855$. *Measured* on a $2000\times2000$ board: density went from $0.2501$ to $0.2858$, in agreement.
 
### 10.3 The mean-field map and why it fails
 
If the independence were preserved after the first generation, the density would follow the iteration $p_{t+1}=F(p_t)$. Its fixed points $F(p)=p$ are
 
| Fixed point | Slope $F'(p)$ | Stability |
|---|---|---|
| $p=0$ | $0$ | stable (extinction) |
| $p\approx0.1925$ | $\approx1.74$ | **unstable** (threshold) |
| $p\approx0.3702$ | $\approx-0.079$ | stable |
 
(A fixed point is stable when $|F'(p)|<1$.) Starting at $0.25$ the mean-field iteration gives $0.2855,\ 0.3294,\ 0.3624,\ 0.3704,\dots\to0.3702$: it predicts the board *fills up* to about 37%.
 
The *measured* density on a $400\times400$ torus does the opposite: $0.2503,\ 0.284,\ 0.2467,\ 0.2337,\ 0.2029$ at generations $0,1,2,5,10$, then $0.118$ (gen. 50), $0.091$ (gen. 100) and $0.063$ (gen. 299), and about $0.05$ after a few hundred generations in the $100\times100$ runs of Section 8.4. From the second step on, neighbouring cells are **correlated** (a cell that just became alive is more likely to sit next to dying cells), so the independence assumption breaks and the mean-field prediction is wrong qualitatively. That is typical of strongly interacting systems, and one of the reasons the long-term density of Life cannot be derived from a closed formula.
 
---
 
## 11. Computability and complexity
 
* **Universality.** The Game of Life is **Turing-complete**: with gliders as signals, guns as clocks and eaters as gates one can build logic circuits and even a universal Turing machine (Berlekamp, Conway and Guy; an explicit machine was built by Rendell). So a Life pattern can compute anything a computer can.
* **Undecidability on the infinite plane.** Because of universality, questions like "will this finite pattern ever die out?" or "will this cell ever become alive?" are **undecidable** on $\mathbb Z^2$: no algorithm can answer them for every starting pattern (reduction from the halting problem).
* **On a torus everything is decidable, in principle.** The state space is finite ($2^{mn}$ states, Theorem 2), so the fate of any starting board can be found by simulating until a repeat. The catch is the length of the orbit, which can be exponential.
* **Cost of one generation in your code.** Each generation costs $\Theta(mn)$ operations (16 `np.roll` calls plus a few element-wise operations; the 4-roll version of Section 3.2 is about 4 times cheaper in rolls but has the same $\Theta(mn)$ order). The computation is a **massively parallel** operation: all cells are updated "at once" with vectorized NumPy.
* **Gosper's Hashlife** uses quadtrees with memoization and can jump exponentially far ahead on regular patterns, at the cost of much more complex code.
---
 
## 12. Other maths hidden in the code
 
**Layout.** For a window of $W\times H$ pixels: $n=\lfloor W/6\rfloor$ columns, $m=\lfloor H/6\rfloor$ rows, and the board is centred with margins $o_x=\lfloor(W-6n)/2\rfloor$ and $o_y=\lfloor(H-6m)/2\rfloor$. Since $0\le W-6n\le5$, the offset is at most $2$ pixels.
 
**Mouse to cell.** A click at pixel $(x,y)$ maps to $c=\lfloor(x-o_x)/6\rfloor$ and $r=\lfloor(y-o_y)/6\rfloor$, and toggling is the map $s\mapsto1-s$, an involution ($1-(1-s)=s$).
 
**Resizing.** `resize_board` copies the overlapping $\min(m,m')\times\min(n,n')$ block and pads the rest with dead cells. The new board is a different torus, so the dynamics change (Section 5).
 
**Time stepping.** With `SPEED = 10` the interval is $\Delta=1000/10=100$ ms. The accumulator $e$ gains $dt$ each frame, and the loop runs $g=\lfloor e/\Delta\rfloor$ generations, then keeps the remainder $e\leftarrow e-g\Delta\in[0,\Delta)$. The average rate is exactly 10 generations per second whatever the frame rate, as long as the computation keeps up. Resetting `elapsed_ms` on pause avoids a burst of generations on resume.
 
**Drawing.** `surfarray` expects an array indexed `[x, y]` while the board is indexed `[row, column]`, hence the transpose `board.T`. A one-pixel-per-cell image is scaled by `CELL_SIZE` with `pygame.transform.scale`, which is a single vectorized resize instead of $mn$ draw calls.
 
---
 
## 13. Pattern catalog
 
`O` is alive, `.` is dead. All patterns below were simulated with your rule; periods and displacements are *measured* values. For the infinite-plane counts (methuselahs) I used a sparse simulation on an unbounded plane, so on your torus the numbers will differ once things wrap around (Section 5).
 
### 13.1 Still lifes (period 1)
 
**Block** (4 cells) | **Beehive** (6) | **Loaf** (7) | **Boat** (5)
 
```
OO     .OO.     .OO.     OO.
OO     O..O     O..O     O.O
       .OO.     .O.O     .O.
                ..O.
```
 
The block is the simplest: each cell has exactly 3 neighbours.
 
### 13.2 Oscillators
 
**Blinker** (3 cells, period 2), the smallest oscillator:
 
```
OOO
```
 
**Toad** (6 cells, period 2):
 
```
.OOO
OOO.
```
 
**Beacon** (8 cells, period 2):
 
```
OO..
OO..
..OO
..OO
```
 
**Pulsar** (48 cells, period 3; it has the full $D_4$ symmetry, so by Section 9.1 it keeps it in every phase):
 
```
..OOO...OOO..
.............
O....O.O....O
O....O.O....O
O....O.O....O
..OOO...OOO..
.............
..OOO...OOO..
O....O.O....O
O....O.O....O
O....O.O....O
.............
..OOO...OOO..
```
 
**Pentadecathlon** (period 15). Start with a row of 10 cells; after about 20 generations it settles into the period-15 oscillator (verified):
 
```
OOOOOOOOOO
```
 
### 13.3 Spaceships
 
**Glider** (5 cells, period 4, moves $(+1,+1)$ per period, speed $c/4$; it moves down-right in the board coordinates of your code):
 
```
.O.
..O
OOO
```
 
**Lightweight spaceship (LWSS)** (9 cells, period 4, moves 2 cells to the left per period, speed $c/2$):
 
```
.O..O
O....
O...O
OOOO.
```
 
### 13.4 Methuselahs: small patterns with a long, rich history
 
All counts below are for the unbounded plane.
 
**R-pentomino** (5 cells). It reaches a maximum of 319 live cells and becomes stable at **generation 1103**, with **116 cells** (including 6 escaped gliders):
 
```
.OO
OO.
.O.
```
 
**Diehard** (8 cells). It lives for a long time and then disappears completely at **generation 130**, with a maximum population of 40:
 
```
......O.
OO......
.O...OOO
```
 
**Acorn** (7 cells). It needs **5206 generations** to stabilize, reaching a maximum of 1057 cells and ending with **633 cells**:
 
```
.O.....
...O...
OO..OOO
```
 
On a small torus these three behave very differently, because their debris wraps around and collides with itself (see the table in Section 5).
 
### 13.5 Gun
 
**Gosper glider gun** (36 cells, period 30). It emits a new glider every 30 generations, so the population grows without bound (linearly), by 5 cells per 30 generations (*measured*: 36, 41, 46, 51 at generations 0, 30, 60, 90). This was the first pattern with unbounded growth (Gosper, 1970):
 
```
........................O...........
......................O.O...........
............OO......OO............OO
...........O...O....OO............OO
OO........O.....O...OO..............
OO........O...O.OO....O.O...........
..........O.....O.......O...........
...........O...O....................
............OO......................
```
 
On a torus, the gliders eventually come around and hit the gun: the infinite growth stops, and the long-term behaviour is whatever the interaction produces.
 
### 13.6 How to place these patterns in your program
 
Your program only lets you click single cells. This helper stamps a pattern on the board, wrapping around the edges consistently with the torus:
 
```python
PATTERNS = {
    "glider":  [".O.", "..O", "OOO"],
    "lwss":    [".O..O", "O....", "O...O", "OOOO."],
    "rpent":   [".OO", "OO.", ".O."],
    "acorn":   [".O.....", "...O...", "OO..OOO"],
    "diehard": ["......O.", "OO......", ".O...OOO"],
}
 
def stamp(board, pattern, row, col):
    """Writes a pattern with its top-left corner at (row, col); wraps like the torus."""
    rows, cols = board.shape
    for i, line in enumerate(pattern):
        for j, ch in enumerate(line):
            if ch == "O":
                board[(row + i) % rows, (col + j) % cols] = 1
```
 
To use it, for example with the `1`-`5` keys while paused, add in `handle_events` (inside the `KEYDOWN` branch):
 
```python
elif self.paused and pygame.K_1 <= event.key <= pygame.K_5:
    name = list(PATTERNS)[event.key - pygame.K_1]
    x, y = pygame.mouse.get_pos()
    col, row = (x - self.offset[0]) // CELL_SIZE, (y - self.offset[1]) // CELL_SIZE
    stamp(self.board, PATTERNS[name], row, col)
```
 
---
 
## 14. Reproducing the experiments
 
All snippets use NumPy only and the same rule as your code.
 
### 14.1 Light cone and perturbation experiment (Section 8)
 
```python
import numpy as np
rng = np.random.default_rng(11)
 
def step(B):                                   # B has shape (trials, m, n)
    nb = np.zeros_like(B)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy or dx:
                nb += np.roll(np.roll(B, dy, 1), dx, 2)
    return ((nb == 3) | ((nb == 2) & (B == 1))).astype(np.int8)
 
T, m, n = 300, 100, 100
A = (rng.random((T, m, n)) < 0.25).astype(np.int8)
P = A.copy()
P[np.arange(T), rng.integers(0, m, T), rng.integers(0, n, T)] ^= 1   # flip one cell per trial
 
x, y = A, P
for t in range(601):
    if t in (0, 1, 5, 10, 25, 50, 100, 200, 400, 600):
        d = (x != y).sum((1, 2))               # Hamming distance per trial
        print(t, f"healed {np.mean(d == 0)*100:.1f}%  median {np.median(d)}")
    x, y = step(x), step(y)
```
 
### 14.2 Exact functional graph of a small torus (Section 6.2)
 
```python
import numpy as np
from collections import Counter
 
def analyze(m, n):
    N = m * n; S = 1 << N
    succ = np.empty(S, dtype=np.int64)
    weights = (1 << np.arange(N, dtype=np.int64))
    CH = 1 << min(N, 20)
    for st in range(0, S, CH):
        ids = np.arange(st, min(S, st + CH), dtype=np.int64)
        B = ((ids[:, None] >> np.arange(N)) & 1).astype(np.int8).reshape(-1, m, n)
        nb = np.zeros_like(B)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dy or dx:
                    nb += np.roll(np.roll(B, dy, 1), dx, 2)
        nxt = ((nb == 3) | ((nb == 2) & (B == 1))).astype(np.int64)
        succ[ids] = (nxt.reshape(-1, N) * weights).sum(1)
 
    indeg = np.bincount(succ, minlength=S)
    gardens_of_eden = int((indeg == 0).sum())
 
    # peel away transient states layer by layer; what remains lies on cycles
    alive = np.ones(S, bool); deg = indeg.copy(); rounds = 0
    frontier = np.where(deg == 0)[0]
    while len(frontier):
        alive[frontier] = False; rounds += 1
        t = succ[frontier]
        deg -= np.bincount(t, minlength=S)
        cand = np.unique(t)
        frontier = cand[(deg[cand] == 0) & alive[cand]]
 
    seen = np.zeros(S, bool); cycles = Counter()
    for s in np.where(alive)[0]:
        if seen[s]: continue
        L, x = 0, s
        while not seen[x]:
            seen[x] = True; x = succ[x]; L += 1
        cycles[L] += 1
    return gardens_of_eden, rounds, dict(sorted(cycles.items()))
 
print(analyze(4, 4))   # (47657, 9, {1: 53, 2: 180, 4: 16, 8: 96})
```
 
(The $5\times5$ case needs about 45 seconds and roughly 1 GB of memory.)
 
### 14.3 Glider period on a torus (Section 9.3)
 
```python
import math
glider = np.array([[0,1,0],[0,0,1],[1,1,1]])
m, n = 20, 30
b = np.zeros((m, n), int); b[1:4, 1:4] = glider
x = b.copy()
for t in range(1, 4 * math.lcm(m, n) + 1):
    x = next_generation(x)
    if np.array_equal(x, b):
        print("returns after", t, "generations"); break     # 240 = 4 * lcm(20, 30)
```
