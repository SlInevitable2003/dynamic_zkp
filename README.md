# dynamic_zkp

这是一个学习项目, 实现 EUROCRYPT 26 上的论文 "Dynamic zk-SNARKs (with Applications to Sparse zk-SNARKs and IVC)".

## 0. 记号与约定

- 固定`m`为2的幂次, `L_i`表示素数域`F`中相应于`m`次单位原根的第`i`个 Lagrange 插值基多项式.
    - 固定`σ`为`(0..m)`上的置换, 则`y_i=L_i-L_{σ^{-1}(i)}`.
    - 承上, `u(X,Y)=Σ{i in 0..m}L_i(X)y_i(Y)`.
- 对于域元素`x`, `[x]_1`和`[x]_2`分别表示`x*g_1`和`x*g_2`, 这里`g_1`和`g_2`分别是群`G_1`和`G_2`的生成元.
- 配对`e: G_1 x G_2 -> G_T`是双线性的.

## 1. Dynamo

### 1.1. Non-Universal Version

#### 1.1.1. Setup

```
func Setup(m, σ: (0..m)上的置换) -> (pk, vk) {
    τ_X, τ_Y <- F;

    for i in 0..m {
        pk.α[i]  = [L_i(τ_X) * (u(τ_X, τ_Y) - y_i(τ_Y)) / (τ_X^m - 1)]_1;
        pk.β[i]  = [y_i(τ_Y) * (L_i(τ_X) - L_i(0)) / τ_X]_1;
        pk.Z[i]  = [τ_Y^m * L_i(τ_X)]_1;
        pk.H[i]  = [τ_X^m * y_i(τ_Y)]_1;
        pk.B[i]  = [y_i(τ_Y) * (L_i(τ_X) - L_i(0)) * τ_X]_1;
        pk.LX[i] = [L_i(τ_X)]_1;
        pk.LY[i] = [L_i(τ_Y)]_1;
    }

    pk.u    = [u(τ_X, τ_Y)]_1;
    pk.X    = [τ_X]_1;
    pk.X2   = [τ_X^2]_1;
    pk.Xm   = [τ_X^m]_1;
    pk.Ym   = [τ_Y^m]_1;
    pk.XmYm = [τ_X^m * τ_Y^m]_1;
    pk.X2Ym = [τ_X^2 * τ_Y^m]_1;

    vk = {
        u    : [u(τ_X, τ_Y)]_2,
        X    : [τ_X]_2,
        X2   : [τ_X^2]_2,
        Xm   : [τ_X^m]_2,
        Ym   : [τ_Y^m]_2,
        invm : [m^{-1}]_2,
    };
    return (pk, vk);
}
```

#### 1.1.2. Prove

```
func Prove(pk, z[0..m], h[0..m]) -> (x, π) {
    ρ_z, ρ_h, ρ_qv, ρ_qh <- F;

    [α]_1 = Σ{i in 0..m} z[i] * pk.α[i];
    [β]_1 = Σ{i in 0..m} z[i] * pk.β[i];
    [Z]_1 = Σ{i in 0..m} z[i] * pk.Z[i];
    [H]_1 = Σ{i in 0..m} z[i] * pk.H[i];
    [B]_1 = Σ{i in 0..m} z[i] * pk.B[i];
    [z]_1 = Σ{i in 0..m} z[i] * pk.LX[i];
    [h]_1 = Σ{i in 0..m} h[i] * pk.LY[i];

    Xm1 = pk.Xm - g_1;
    Ym1 = pk.Ym - g_1;

    [z_zk]_1 = [z]_1 + ρ_z * Xm1;
    [h_zk]_1 = [h]_1 + ρ_h * Ym1;

    [q_zk]_1  = (-ρ_h / m) * g_1 + ρ_qv * Xm1 + ρ_qh * pk.X;
    [α_zk]_1  = [α]_1 + ρ_z * pk.u - ρ_qv * Ym1;
    [β_zk]_1  = [β]_1 - ρ_qh * Ym1;
    [Z_zk]_1  = [Z]_1 + ρ_z * (pk.XmYm - pk.Ym);
    [H_zk]_1  = [H]_1 + ρ_h * (pk.XmYm - pk.Xm);
    [B_zk]_1  = [B]_1 - ρ_qh * (pk.X2Ym - pk.X2);

    x = ([z_zk]_1, [h_zk]_1);
    π = ([α_zk]_1, [β_zk]_1, [q_zk]_1, [Z_zk]_1, [H_zk]_1, [B_zk]_1);
    return (x, π);
}
```

#### 1.1.3. Verify

```
func Verify(vk, x, π) -> bool {
    ([z_zk]_1, [h_zk]_1) = x;
    ([α_zk]_1, [β_zk]_1, [q_zk]_1, [Z_zk]_1, [H_zk]_1, [B_zk]_1) = π;

    Xm1 = vk.Xm - g_2;
    Ym1 = vk.Ym - g_2;

    if e([z_zk]_1, vk.u) != e([α_zk]_1, Xm1) + e([β_zk]_1, vk.X) + e([h_zk]_1, vk.invm) + e([q_zk]_1, Ym1) { return false; }
    if e([z_zk]_1, vk.Ym) != e([Z_zk]_1, g_2) { return false; }
    if e([h_zk]_1, vk.Xm) != e([H_zk]_1, g_2) { return false; }
    if e([β_zk]_1, vk.X2) != e([B_zk]_1, g_2) { return false; }
    return true;
}
```