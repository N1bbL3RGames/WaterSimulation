using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FluidSimluation : MonoBehaviour
{
    enum Cell { EMPTY, FLUID, SOLID };
    const int N = 32;
    const int Np1 = N + 1;
    const int Y = N * (N + 1);
    const int PN = (N - 2) / 4 * (N - 4) * 4;

    double[] u = new double[2 * N * (N + 1)];
    double[] ux = new double[2 * N * (N + 1)];

    double[] m = new double[2 * N * (N + 1)];
    double[] p = new double[N * N];

    double[] divu = new double[N * N];

    Cell[] flags = new Cell[N * N];

    double[] px = new double[PN], py = new double[PN], pvx = new double[PN], pvy = new double[PN];

    //double dt = 0.01; //Time.deltaTime;
    public GameObject particle;
    public GameObject wall;
    GameObject[] particles = new GameObject[PN];
    GameObject[] walls = new GameObject[4 * (N - 1)];

    // Start is called before the first frame update
    void Start()
    {
        int idx = 0;
        for (int a = 0; a < pvx.Length; a++)
        {
            pvx[a] = 0;
            pvy[a] = 0;
        }

        for (int j = 1; j < N - 1; ++j)
            for (int i = 1; i < N - 1; ++i)
                if (idx < PN && i - 1 < (N - 2) / 4)
                    for (int jj = 0; jj < 2; ++jj)
                        for (int ii = 0; ii < 2; ++ii)
                        {
                            px[idx] = i + 0.25 + ii * 0.5;
                            py[idx++] = j + 0.25 + jj * 0.5;
                        }

        for (int a = 0; a < PN; a++)
            particles[a] = Instantiate(particle, new Vector3((float)px[a], (float)py[a], 0), Quaternion.identity);
    }

    void particles2grid()
    {
        for (int a = 0; a < u.Length; a++)
        {
            u[a] = 0;
            m[a] = 0;
        }

        for (int j = 0; j < N; ++j) 
            for (int i = 0; i < N; ++i)
                flags[i + j * N] = (i == 0 || j == 0 || i == N - 1 || j == N - 1) ? Cell.SOLID : Cell.EMPTY;

        for (int k = 0; k < PN; ++k)
        {
            int i = (int)px[k], j = (int)py[k], fi = (int)(px[k] - 0.5), fj = (int)(py[k] - 0.5);
            flags[i + j * N] = Cell.FLUID;

            for (int jj = fj; jj < fj + 2; ++jj) 
                for (int ii = i; ii < i + 2; ++ii)
                {
                    u[ii + jj * Np1] += pvx[k] * (1 - Mathf.Abs((float)(ii - px[k]))) * (1 - Mathf.Abs((float)(jj + 0.5 - py[k])));
                    m[ii + jj * Np1] += (1 - Mathf.Abs((float)(ii - px[k]))) * (1 - Mathf.Abs((float)(jj + 0.5 - py[k])));
                }

            for (int jj = j; jj < j + 2; ++jj) 
                for (int ii = fi; ii < fi + 2; ++ii)
                {
                    u[Y + ii + jj * N] += pvy[k] * (1 - Mathf.Abs((float)(ii + 0.5 - px[k]))) * (1 - Mathf.Abs((float)(jj - py[k])));
                    m[Y + ii + jj * N] += (1 - Mathf.Abs((float)(ii + 0.5 - px[k]))) * (1 - Mathf.Abs((float)(jj - py[k])));
                }
        }

        for (int k = 0; k < 2 * N * Np1; ++k)
            if (m[k] > 1e-8) 
            { 
                m[k] = 1 / m[k]; 
                u[k] *= m[k]; 
            }

        ux = u;
    }

    void pressureSolve()
    {
        for (int j = 1; j < N - 1; ++j) 
            for (int i = 1; i < N - 1; ++i)
                divu[i + j * N] = -(u[i + 1 + j * Np1] - u[i + j * Np1] + u[Y + i + (j + 1) * N] - u[Y + i + j * N]);

        for (int k = 1; k < N - 1; ++k)
        {
            divu[1 + k * N] -= u[1 + k * Np1];
            divu[N - 2 + k * N] += u[N - 1 + k * Np1];
            divu[k + N] -= u[Y + k + N];
            divu[k + (N - 2) * N] += u[Y + k + (N - 1) * N];
        }

        for (int a = 0; a < p.Length; a++)
            p[a] = 0;

        for (int k = 0; k < 10 * N * N; ++k)    // Gauss-Seidel solve for pressure
            for (int j = 1; j < N - 1; ++j) 
                for (int i = 1; i < N - 1; ++i)
                {
                    if (flags[i + j * N] != Cell.FLUID) 
                        continue;

                    double l = i > 1 ? -1 : 0, r = i < N - 2 ? -1 : 0, b = j > 1 ? -1 : 0, t = j < N - 2 ? -1 : 0, c = -(l + r + b + t);
                    p[i + j * N] = 1 / c * (divu[i + j * N] - l * p[i - 1 + j * N] - r * p[i + 1 + j * N] - b * p[i + j * N - N] - t * p[i + j * N + N]);
                }

        for (int j = 1; j < N - 1; ++j)
            for (int i = 1; i < N; ++i)
                u[i + j * Np1] -= p[i + j * N] - p[i - 1 + j * N];

        for (int j = 1; j < N; ++j)
            for (int i = 1; i < N - 1; ++i)
                u[Y + i + j * N] -= p[i + j * N] - p[i + j * N - N];

        for (int k = 1; k < N - 1; ++k)
            u[1 + k * Np1] = u[N - 1 + k * Np1] = u[Y + k + N] = u[Y + k + (N - 1) * N] = 0;
    }

    void updateAndAdvect(double flip)
    {
        for (int k = 0; k < PN; ++k)
        {
            int i = (int)px[k], j = (int)py[k], fi = (int)(px[k] - 0.5), fj = (int)(py[k] - 0.5);
            double vx = 0, vy = 0, vxOld = 0, vyOld = 0;

            for (int jj = fj; jj < fj + 2; ++jj) 
                for (int ii = i; ii < i + 2; ++ii)
                {
                    vx += u[ii + jj * Np1] * (1 - Mathf.Abs((float)(ii - px[k]))) * (1 - Mathf.Abs((float)(jj + 0.5 - py[k])));
                    vxOld += ux[ii + jj * Np1] * (1 - Mathf.Abs((float)(ii - px[k]))) * (1 - Mathf.Abs((float)(jj + 0.5 - py[k])));
                }

            for (int jj = j; jj < j + 2; ++jj) 
                for (int ii = fi; ii < fi + 2; ++ii)
                {
                    vy += u[Y + ii + jj * N] * (1 - Mathf.Abs((float)(ii + 0.5 - px[k]))) * (1 - Mathf.Abs((float)(jj - py[k])));
                    vyOld += ux[Y + ii + jj * N] * (1 - Mathf.Abs((float)(ii + 0.5 - px[k]))) * (1 - Mathf.Abs((float)(jj - py[k])));
                }

            pvx[k] = (1 - flip) * vx + flip * (pvx[k] + vx - vxOld);
            pvy[k] = (1 - flip) * vy + flip * (pvy[k] + vy - vyOld);
            px[k] = (double)Mathf.Min(Mathf.Max((float)(px[k] + vx * Time.deltaTime), (float)1.001), (float)(N - 1.001));
            py[k] = (double)Mathf.Min(Mathf.Max((float)(py[k] + vy * Time.deltaTime), (float)1.001), (float)(N - 1.001));
        }
    }

    // Update is called once per frame
    void Update()
    {
        particles2grid();

        for (int k = 0; k < N * Np1; ++k)
            if (m[Y + k] > 1e-8) 
                u[Y + k] += -9.81 * Time.deltaTime * (N - 2);

        pressureSolve();
        updateAndAdvect(0.96);

        float widthScalar = 0.16f;
        float heightScalar = 0.9f;

        for (int k = 0; k < PN; ++k)
            particles[k].gameObject.transform.position = new Vector3((float)px[k] * widthScalar * 16 - 8, (float)py[k] * heightScalar * 10 - 5, 0);
        /*
        int idx = 0;
        for (int k = 0; k < walls.Length; k++)
            Destroy(walls[k]);

        for (int j = 0; j < N; ++j) 
            for (int i = 0; i < N; ++i)
            {
                if (flags[i + j * N] != Cell.SOLID) 
                    continue;

                walls[idx++] = Instantiate(wall, new Vector3((float)(1.9 * (i / (double)N - 0.5)), (float)(1.9 * (j / (double)N - 0.5)), 0), Quaternion.identity);
                walls[idx++] = Instantiate(wall, new Vector3((float)(1.9 * (i / (double)N - 0.5)), (float)(1.9 * ((j + 1) / (double)N - 0.5)), 0), Quaternion.identity);
                walls[idx++] = Instantiate(wall, new Vector3((float)(1.9 * ((i + 1) / (double)N - 0.5)), (float)(1.9 * ((j + 1) / (double)N - 0.5)), 0), Quaternion.identity);
                walls[idx++] = Instantiate(wall, new Vector3((float)(1.9 * ((i + 1) / (double)N - 0.5)), (float)(1.9 * (j / (double)N - 0.5)), 0), Quaternion.identity);
            }*/
    }
}
