void LaxWendroffStep() {
    // compute flux F from U
    for (int j = 0; j < N; j++) {
        double rho = U[j][0];
        double m = U[j][1];
        double e = U[j][2];
        double p = (gama - 1) * (e - m * m / rho / 2);
        F[j][0] = m;
        F[j][1] = m * m / rho + p;
        F[j][2] = m / rho * (e + p);
    }
    // half step
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            newU[j][i] = (U[j + 1][i] + U[j][i]) / 2 -
                         tau / 2 / h * (F[j + 1][i] - F[j][i]);
    boundaryConditions(newU);

    // compute flux at half steps
    for (int j = 0; j < N; j++) {
        double rho = newU[j][0];
        double m = newU[j][1];
        double e = newU[j][2];
        double p = (gama - 1) * (e - m * m / rho / 2);
        F[j][0] = m;
        F[j][1] = m * m / rho + p;
        F[j][2] = m / rho * (e + p);
    }
    // step using half step flux
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            newU[j][i] = U[j][i] - tau / h * (F[j][i] - F[j - 1][i]);
    // update U from newU
    for (int j = 1; j < N - 1; j++)
        for (int i = 0; i < 3; i++)
            U[j][i] = newU[j][i];
}
