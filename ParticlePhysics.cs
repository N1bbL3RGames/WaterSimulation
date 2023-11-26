using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public class ParticlePhysics : MonoBehaviour
{
    // Scene variables
    float gravity = 9.81f;
    float flipRatio = 0.9f;
    int numPressureIterations = 50;
    int numParticleIterations = 2;
    float relaxivity = 1.9f;
    bool compensateDrift = true;
    bool separateParticles = true;
    Vector2 objectPosition = Vector2.zero;
    float objectRadius = 0.15f;
    bool paused = false;
    Vector2 objectVelocity = Vector2.zero;
    bool showParticles = true;
    FlipFluid fluid = null;
    //int frameNumber = 0;
    //bool showObstacle = true;
    //bool showGrid = false;

    // Unity variables
    static float tankWidth = 4;
    static float tankHeight = 3;

    public float density = 1000f;
    public float spacing = UnityEngine.Screen.height / 100f;
    public float particleRadius = 0.3f * UnityEngine.Screen.height / 100f;
    public GameObject particlePrefab;

    GameObject[] particles;

    void Start()
    {
        SetupScene();
    }

    void Update()
    {
        Simulate();
        Draw();
    }

    void SetupScene()
    {
        float dx = 2.0f * particleRadius;
        float dy = Mathf.Sqrt(3.0f) / 2.0f * dx;
        float x = Mathf.Floor((0.6f * UnityEngine.Screen.width - 2.0f * spacing - 2.0f * particleRadius) / dx);
        float y = Mathf.Floor((0.8f * UnityEngine.Screen.height - 2.0f * spacing - 2.0f * particleRadius) / dy);
        int maxParticles = (int)(x * y);
        Debug.Log(maxParticles);

        fluid = new FlipFluid(density, spacing, particleRadius, maxParticles);
        fluid.numParticles = (int)(x * y);
        int index = 0;

        for (int a = 0; a < x; a++)
            for (int b = 0; b < y; b++)
                fluid.particlePosition[index++] = new Vector2(spacing + particleRadius + dx * a + (b % 2 == 0 ? 0 : particleRadius), spacing + particleRadius + dy * b);

        for (int a = 0; a < fluid.flip.x; a++)
            for (int b = 0; b < fluid.flip.y; b++)
                fluid.s[(int)(a * fluid.flip.y + b)] = (a == 0 || a == fluid.flip.x - 1 || b == 0) ? 0.0f : 1.0f;

        SetObject(new Vector2(3.0f, 2.0f), true);
        particles = new GameObject[maxParticles];

        for (int a = 0; a < fluid.numParticles; a++)
            particles[a] = Instantiate(particlePrefab, new Vector3(fluid.particlePosition[a].x / UnityEngine.Screen.width * 16 - 8, fluid.particlePosition[a].y / UnityEngine.Screen.height * 10 - 5, 0), Quaternion.identity);
    }

    Vector2 testing = Vector2.zero;

    private void Draw()
    {
        if (showParticles)
        {
            float size = 2.0f * fluid.particleRadius;

            for (int a = 0; a < particles.Length; a++) // X: -8 to 8, Y: -5 to 5
                particles[a].gameObject.transform.position = new Vector3(fluid.particlePosition[a].x / UnityEngine.Screen.width * 16 - 8, fluid.particlePosition[a].y / UnityEngine.Screen.height * 10 - 5, 0);
        }
    }

    void SetObject(Vector2 position, bool reset)
    {
        Vector2 velocity = Vector2.zero;

        if (!reset)
            velocity = new Vector2((position.x - objectPosition.x) / Time.deltaTime, (position.y - objectPosition.y) / Time.deltaTime);

        objectPosition = position;

        float cd = Mathf.Sqrt(2) * fluid.cellHeight;

        for (int a = 1; a < fluid.flip.x - 2; a++)
            for (int b = 1; b < fluid.flip.y - 2; b++)
            {
                fluid.s[(int)(a * fluid.flip.y + b)] = 1.0f;
                float dx = (a + 0.5f) * fluid.cellHeight - position.x;
                float dy = (b + 0.5f) * fluid.cellHeight - position.y;

                if (dx * dx + dy * dy < objectRadius * objectRadius)
                {
                    fluid.s[(int)(a * fluid.flip.y + b)] = 0;
                    fluid.u[(int)(a * fluid.flip.y + b)] = velocity.x;
                    fluid.u[(int)((a + 1) * fluid.flip.y + b)] = velocity.x;
                    fluid.v[(int)(a * fluid.flip.y + b)] = velocity.y;
                    fluid.v[(int)((a + 1) * fluid.flip.y + b)] = velocity.y;
                }
            }

        //showObstacle = true;
        objectVelocity = velocity;
    }

    void Simulate()
    {
        if (!paused)
        {
            fluid.Simulate(gravity, flipRatio, numPressureIterations, numParticleIterations, relaxivity, compensateDrift, separateParticles, objectPosition, objectVelocity, objectRadius);
            //frameNumber++;
        }
    }

    class FlipFluid
    {
        // Constructor variables
        public float density;
        public float spacing;
        public int maxParticles;
        public float particleRadius;

        // Simulation variables
        public Vector2 flip;
        public float cellHeight, borderSpacing;
        public int numCells;

        public float[] u, v, du, dv, oldU, oldV, p, s;
        public CellType[] cell;
        public Color[] cellColor;

        public Vector2[] particlePosition, particleVelocity;
        public Color[] particleColor;
        public float[] particleDensity;
        public float particleRestingDensity;
        public float particleSpacing;
        public Vector2 particle;
        public int particleCells;

        public int[] numCellParticle;
        public int[] firstCellParticle;
        public int[] cellParticleID;

        public int numParticles;

        public enum CellType
        {
            Empty, Fluid, Solid
        };

        public FlipFluid(float density, float spacing, float particleRadius, int maxParticles)
        {
            this.density = density;
            this.spacing = spacing;
            this.maxParticles = maxParticles;
            this.particleRadius = particleRadius;

            this.flip = new Vector2(Mathf.Floor(UnityEngine.Screen.width / spacing) + 1, Mathf.Floor(UnityEngine.Screen.height / spacing) + 1);
            this.cellHeight = Mathf.Max(UnityEngine.Screen.width / flip.x, UnityEngine.Screen.height / flip.y);
            this.borderSpacing = 1.0f / cellHeight;
            this.numCells = (int)flip.x * (int)flip.y;

            this.u = new float[numCells];
            this.v = new float[numCells];
            this.du = new float[numCells];
            this.dv = new float[numCells];
            this.oldU = new float[numCells];
            this.oldV = new float[numCells];
            this.p = new float[numCells];
            this.s = new float[numCells];

            this.cell = new CellType[numCells];
            this.cellColor = new Color[numCells];
            this.particlePosition = new Vector2[maxParticles];
            this.particleColor = new Color[maxParticles];

            for (int a = 0; a < particleColor.Length; a++)
                particleColor[a] = Color.blue;

            this.particleVelocity = new Vector2[maxParticles];
            this.particleDensity = new float[numCells];
            this.particleRestingDensity = 0;
            this.particleSpacing = 1.0f / (2.2f * particleRadius);
            this.particle = new Vector2(Mathf.Floor(UnityEngine.Screen.width * particleSpacing) + 1, Mathf.Floor(UnityEngine.Screen.height * particleSpacing) + 1);
            this.particleCells = (int)particle.x * (int)particle.y;

            this.numCellParticle = new int[particleCells];
            this.firstCellParticle = new int[particleCells + 1];
            this.cellParticleID = new int[maxParticles];

            this.numParticles = 0;
        }

        void IntegrateParticles(float gravity)
        {
            for (int a = 0; a < numParticles; a++)
            {
                particleVelocity[a] += Vector2.down * Time.deltaTime * gravity;
                particlePosition[a] += particleVelocity[a] * Time.deltaTime;
            }
        }

        void PushParticlesApart(int numIterations)
        {
            float diffusionCoefficient = 0.001f;

            for (int a = 0; a < numCellParticle.Length; a++)
                numCellParticle[a] = 0;

            for (int a = 0; a < numParticles; a++)
            {
                float x = Mathf.Clamp(Mathf.Floor(particlePosition[a].x * particleSpacing), 0, particle.x - 1);
                float y = Mathf.Clamp(Mathf.Floor(particlePosition[a].y * particleSpacing), 0, particle.y - 1);
                numCellParticle[(int)(x * particle.y + y)]++;
            }

            int firstNum = 0;

            for (int a = 0; a < particleCells; a++)
            {
                firstNum += numCellParticle[a];
                firstCellParticle[a] = firstNum;
            }

            firstCellParticle[particleCells] = firstNum;

            for (int a = 0; a < numParticles; a++)
            {
                float x = Mathf.Clamp(Mathf.Floor(particlePosition[a].x * particleSpacing), 0, particle.x - 1);
                float y = Mathf.Clamp(Mathf.Floor(particlePosition[a].y * particleSpacing), 0, particle.y - 1);
                firstCellParticle[(int)(x * particle.y + y)]--;
                cellParticleID[firstCellParticle[(int)(x * particle.y + y)]] = a;
            }

            float forceRadius = 2.0f * particleRadius;
            float forceRadiusSqr = Mathf.Pow(forceRadius, 2);

            for (int a = 0; a < numIterations; a++)
                for (int b = 0; b < numParticles; b++)
                {
                    float x = Mathf.Floor(particlePosition[b].x * particleSpacing);
                    float y = Mathf.Floor(particlePosition[b].y * particleSpacing);

                    Vector2 v0 = new Vector2(Mathf.Max(x - 1, 0), Mathf.Max(y - 1, 0));
                    Vector2 v1 = new Vector2(Mathf.Max(x + 1, particle.x - 1), Mathf.Max(y - 1, particle.y - 1));

                    for (int c = (int)v0.x; c <= (int)v1.x; c++)
                        for (int d = (int)v0.y; d <= (int)v1.y; d++)
                        {
                            int cellNumber = (int)(c * particle.y + d);
                            int first = firstCellParticle[cellNumber % firstCellParticle.Length];
                            int last = firstCellParticle[(cellNumber + 1) % firstCellParticle.Length];

                            for (int e = first; e < last; e++)
                            {
                                int id = cellParticleID[e];

                                if (id == b)
                                    continue;

                                float dx = particlePosition[id].x - particlePosition[b].x;
                                float dy = particlePosition[id].y - particlePosition[b].y;
                                float dSqr = dx * dx + dy * dy;

                                if (dSqr > forceRadiusSqr || dSqr == 0)
                                    continue;

                                float derivative = Mathf.Sqrt(dSqr);
                                float s = (forceRadius - derivative) / (2 * derivative);
                                dx *= s;
                                dy *= s;

                                particlePosition[b] -= new Vector2(dx, dy);
                                particlePosition[id] += new Vector2(dx, dy);

                                Color diffusedColor = (particleColor[b] + particleColor[id]) / 2;
                                particleColor[b] = particleColor[b] + (diffusedColor - particleColor[b]) * diffusionCoefficient;
                                particleColor[id] = particleColor[id] + (diffusedColor - particleColor[id]) * diffusionCoefficient;
                            }
                        }
                }
        }

        void HandleParticleCollisions(Vector2 objectPosition, Vector2 objectVelocity, float objectRadius)
        {
            float h = 1.0f / borderSpacing;
            float orSqr = objectRadius * objectRadius;

            float forceRadius = objectRadius + particleRadius;
            float forceRadiusSqr = Mathf.Pow(forceRadius, 2);

            float minXY = h + particleRadius;
            float maxX = (flip.x - 1) * h - particleRadius;
            float maxY = (flip.y - 1) * h - particleRadius;

            for (int a = 0; a < numParticles; a++)
            {
                float dx = particlePosition[a].x - objectPosition.x;
                float dy = particlePosition[a].y - objectPosition.y;
                float dSqr = dx * dx + dy * dy;

                if (dSqr < forceRadiusSqr)
                    particleVelocity[a] = objectVelocity;

                if (particlePosition[a].x < minXY)
                {
                    particlePosition[a].x = minXY;
                    particleVelocity[a].x = 0;
                }

                if (particlePosition[a].x > maxX)
                {
                    particlePosition[a].x = maxX;
                    particleVelocity[a].x = 0;
                }

                if (particlePosition[a].y < minXY)
                {
                    particlePosition[a].y = minXY;
                    particleVelocity[a].y = 0;
                }

                if (particlePosition[a].y > maxY)
                {
                    particlePosition[a].y = maxY;
                    particleVelocity[a].y = 0;
                }
            }
        }

        void UpdateParticleDensity()
        {
            float halfBorder = borderSpacing / 2;

            for (int a = 0; a < particleDensity.Length; a++)
                particleDensity[a] = 0;

            for (int a = 0; a < numParticles; a++)
            {
                float x = Mathf.Clamp(particlePosition[a].x, cellHeight, (flip.x - 1) * cellHeight);
                float y = Mathf.Clamp(particlePosition[a].y, cellHeight, (flip.y - 1) * cellHeight);

                Vector2 v0 = new Vector2(Mathf.Floor((particlePosition[a].x - halfBorder) * borderSpacing), Mathf.Floor((particlePosition[a].y - halfBorder) * borderSpacing));
                Vector2 v1 = new Vector2(Mathf.Min(v0.x + 1, flip.x - 2), Mathf.Min(v0.y + 1, flip.y - 2));

                Vector2 vt = new Vector2(((particlePosition[a].x - halfBorder) - v0.x * cellHeight) * borderSpacing, ((particlePosition[a].y - halfBorder) - v0.y * cellHeight) * borderSpacing);
                Vector2 vs = new Vector2(1.0f - vt.x, 1.0f - vt.y);

                if (v0.x < flip.x && v0.y < flip.y)
                    particleDensity[(int)(v0.x * flip.y + v0.y)] += vs.x * vs.y;

                if (v1.x < flip.x && v0.y < flip.y)
                    particleDensity[(int)(v1.x * flip.y + v0.y)] += vt.x * vs.y;

                if (v1.x < flip.x && v1.y < flip.y)
                    particleDensity[(int)(v1.x * flip.y + v1.y)] += vt.x * vt.y;

                if (v0.x < flip.x && v1.y < flip.y)
                    particleDensity[(int)(v0.x * flip.y + v1.y)] += vs.x * vt.y;
            }

            if (particleRestingDensity == 0)
            {
                float sum = 0;
                float numCells = 0;

                for (int a = 0; a < numCells; a++)
                    if (cell[a] == CellType.Fluid)
                    {
                        sum += particleDensity[a];
                        numCells++;
                    }

                if (numCells > 0)
                    particleRestingDensity = sum / numCells;
            }
        }

        void TransferVelocities(bool toCell, float flipRatio)
        {
            float halfBorder = borderSpacing / 2;

            if (toCell)
            {
                u.CopyTo(oldU, 0);
                v.CopyTo(oldV, 0);

                for (int a = 0; a < du.Length; a++)
                {
                    du[a] = 0;
                    dv[a] = 0;
                    u[a] = 0;
                    v[a] = 0;
                }

                for (int a = 0; a < numCells; a++)
                    cell[a] = s[a] == 0 ? CellType.Solid : CellType.Empty;

                for (int a = 0; a < numParticles; a++)
                {
                    float x = Mathf.Clamp(Mathf.Floor(particlePosition[a].x * borderSpacing), 0, flip.x - 1);
                    float y = Mathf.Clamp(Mathf.Floor(particlePosition[a].y * borderSpacing), 0, flip.y - 1);
                    if (cell[(int)(x * flip.y + y)] == CellType.Empty)
                        cell[(int)(x * flip.y + y)] = CellType.Fluid;
                }
            }
            
            for (int a = 0; a < 2; a++)
            {
                float dx = a == 0 ? 0 : halfBorder;
                float dy = a == 0 ? halfBorder : 0;

                float[] f = a == 0 ? u : v;
                float[] oldF = a == 0 ? oldU : oldV;
                float[] d = a == 0 ? du : dv;

                for (int b = 0; b < particleVelocity.Length; b++)
                {
                    float x = Mathf.Clamp(particlePosition[b].x, cellHeight, (flip.x - 1) * cellHeight);
                    float y = Mathf.Clamp(particlePosition[b].y, cellHeight, (flip.y - 1) * cellHeight);

                    Vector2 v0 = new Vector2(Mathf.Min(Mathf.Floor((x - dx) * borderSpacing), flip.x - 2), Mathf.Min(Mathf.Floor((y - dy) * borderSpacing), flip.y - 2));
                    Vector2 v1 = new Vector2(Mathf.Min(v0.x + 1, flip.x - 2), Mathf.Min(v0.y + 1, flip.y - 2));

                    Vector2 vt = new Vector2(((x - dx) - v0.x * cellHeight) * borderSpacing, ((y - dy) - v0.y * cellHeight) * borderSpacing);
                    Vector2 vs = new Vector2(1.0f - vt.x, 1.0f - vt.y);

                    float[] dr = { vs.x * vs.y, vt.x * vs.y, vt.x * vt.y, vs.x * vt.y };
                    float[] nr = { v0.x * flip.y + v0.y, v1.x * flip.y + v0.y, v1.x * flip.y + v1.y, v0.x * flip.y + v1.y };
                    //Debug.Log(v0.x + " " + v0.y + " " + v1.x + " " + v1.y + " " + flip.y);
                    Debug.Log(nr[0] + " " + nr[1] + " " + nr[2] + " " + nr[3]);

                    float velocity = a == 0 ? particleVelocity[b].x : particleVelocity[b].y;
                    
                    if (toCell)
                        for (int c = 0; c <= 3; c++)
                        {
                            f[(int)nr[c]] += velocity * dr[c];
                            d[(int)nr[c]] += dr[c];
                            //try
                            //{
                                
                            //}
                            //catch(Exception e)
                            //{
                            //    Debug.Log(f.Length + " : " + (int)nr[c]);
                            //}
                        }
                    else
                    {
                        float offset = a == 0 ? flip.y : 1;
                        float[] valid = new float[4];

                        for (int c = 0; c < valid.Length; c++)
                            valid[c] = cell[(int)nr[c]] != CellType.Empty || cell[(int)(nr[c] - offset)] != CellType.Empty ? 1.0f : 0.0f;

                        float d2 = valid[0] * dr[0] + valid[1] * dr[1] + valid[2] * dr[2] + valid[3] * dr[3];

                        if (d2 > 0.0f)
                        {
                            float pic = (valid[0] * dr[0] * f[(int)nr[0]] + valid[1] * dr[1] * f[(int)nr[1]] + valid[2] * dr[2] * f[(int)nr[2]] + valid[3] * dr[3] * f[(int)nr[3]]) / d2;
                            float cor = (valid[0] * dr[0] * (f[(int)nr[0]] - oldF[(int)nr[0]]) + valid[1] * dr[1] * (f[(int)nr[1]] - oldF[(int)nr[1]]) + valid[2] * dr[2] * (f[(int)nr[2]] - oldF[(int)nr[2]]) + valid[3] * dr[3] * (f[(int)nr[3]] - oldF[(int)nr[3]])) / d2;
                            float flipV = velocity * cor;

                            if (a == 0)
                                particleVelocity[b].x = (1.0f - flipRatio) * pic + flipRatio * flipV;
                            else
                                particleVelocity[b].y = (1.0f - flipRatio) * pic + flipRatio * flipV;
                        }
                    }
                }

                if (toCell)
                {
                    for (int b = 0; b < f.Length; b++)
                        if (d[b] > 0)
                            f[b] /= d[b];

                    for (int b = 0; b < flip.x; b++)
                        for (int c = 0; c < flip.y; c++)
                        {
                            bool solid = cell[(int)(b * flip.y + c)] == CellType.Solid;

                            if (solid || (b > 0 && cell[(int)((b - 1) * flip.y + c)] == CellType.Solid))
                                u[(int)(b * flip.y + c)] = oldU[(int)(b * flip.y + c)];
                            if (solid || (c > 0 && cell[(int)(b * flip.y + c - 1)] == CellType.Solid))
                                v[(int)(b * flip.y + c)] = oldV[(int)(b * flip.y + c)];
                        }
                }
            }
        }

        void SolveIncompressibility(int numIterations, float relaxivity, bool compensateDrift)
        {
            for (int a = 0; a < p.Length; a++)
                p[a] = 0;

            u.CopyTo(oldU, 0);
            v.CopyTo(oldV, 0);

            float cp = density * cellHeight / Time.deltaTime;

            for (int a = 0; a < numIterations; a++)
                for (int b = 1; b < flip.x - 1; b++)
                    for (int c = 1; c < flip.y - 1; c++)
                    {
                        if (cell[(int)(b * flip.y + c)] != CellType.Fluid)
                            continue;

                        int center = (int)(b * flip.y + c);
                        int left = (int)((b - 1) * flip.y + c);
                        int right = (int)((b + 1) * flip.y + c);
                        int bottom = (int)(b * flip.y + c - 1);
                        int top = (int)(b * flip.y + c + 1);

                        float ss = s[left] + s[right] + s[bottom] + s[top];

                        if (ss == 0)
                            continue;

                        float div = u[right] - u[center] + v[top] - v[center];

                        if (particleRestingDensity > 0 && compensateDrift)
                        {
                            float compression = particleDensity[(int)(b * flip.y + c)] - particleRestingDensity;

                            if (compression > 0)
                                div -= compression;
                        }

                        float rho = -div / ss * relaxivity;
                        p[center] += cp * rho;

                        u[center] -= s[left] * rho;
                        u[right] += s[right] * rho;
                        v[center] -= s[bottom] * rho;
                        v[top] += s[top] * rho;
                    }
        }

        void UpdateParticleColors()
        {
            for (int a = 0; a < numParticles; a++)
            {
                float step = 0.01f;

                particleColor[a] = new Color(Mathf.Clamp(particleColor[a].r - step, 0, 1), Mathf.Clamp(particleColor[a].g - step, 0, 1), Mathf.Clamp(particleColor[a].b + step, 0, 1));

                float x = Mathf.Clamp(Mathf.Floor(particlePosition[a].x * borderSpacing), 1, flip.x - 1);
                float y = Mathf.Clamp(Mathf.Floor(particlePosition[a].y * borderSpacing), 1, flip.y - 1);

                if (particleRestingDensity > 0)
                {
                    float relativeDensity = particleDensity[(int)(x * flip.y + y)] / particleRestingDensity;

                    if (relativeDensity < 0.7)
                        particleColor[a] = new Color(0.8f, 0.8f, 1.0f);
                }
            }
        }

        void SetSciColor(int cellNumber, float value, float min, float max)
        {
            value = (max - min) == 0 ? 0.5f : (Mathf.Min(Mathf.Max(value, min), max - 0.0001f) - min) / (max - min);
            float num = Mathf.Floor(value * 4);
            float sci = value * 4 + -num;

            switch (num)
            {
                case 0:
                    cellColor[cellNumber] = new Color(0, sci, 1);
                    break;
                case 1:
                    cellColor[cellNumber] = new Color(0, 1, 1 - sci);
                    break;
                case 2:
                    cellColor[cellNumber] = new Color(sci, 1, 0);
                    break;
                case 3:
                    cellColor[cellNumber] = new Color(1, 1 - sci, 0);
                    break;
            }
        }

        void UpdateCellColors()
        {
            for (int a = 0; a < cellColor.Length; a++)
                cellColor[a] = new Color(0, 0, 0);

            for (int a = 0; a < numCells; a++)
            {
                if (cell[a] == CellType.Solid)
                    cellColor[a] = new Color(0.5f, 0.5f, 0.5f);
                else if (cell[a] == CellType.Fluid)
                {
                    float rho = particleDensity[a];

                    if (particleRestingDensity > 0)
                        rho /= particleRestingDensity;

                    SetSciColor(a, rho, 0, 2);
                }
            }
        }

        public void Simulate(float gravity, float flipRatio, int numPressureIterations, int numParticleIterations, float relaxivity, bool compensateDrift, bool separateParticles, Vector2 objectPosition, Vector2 objectVelocity, float objectRadius)
        {
            for (int a = 0; a < Time.deltaTime; a++)
            {
                IntegrateParticles(gravity);

                if (separateParticles)
                    PushParticlesApart(numParticleIterations);
                
                HandleParticleCollisions(objectPosition, objectVelocity, objectRadius);
                TransferVelocities(true, flipRatio);
                UpdateParticleDensity();
                SolveIncompressibility(numPressureIterations, relaxivity, compensateDrift);
                TransferVelocities(false, flipRatio);
            }

            //UpdateParticleColors();
            //UpdateCellColors();
        }
    }
}