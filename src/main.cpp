
#include <mpi.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// #include <GL/glew.h>
// #include <GLFW/glfw3.h>
// #include <glm/glm.hpp>
// #include <GL/glu.h>
// #include <GL/glut.h>

const double G = 0.0001;
const double rlimit = 0.03;

const int MAX_X = 4;
const int MAX_Y = 4;

struct options_t
{
  const char *in_file;
  const char *out_file;
  int steps;
  double theta;
  double timestep;
  bool visualization;
};
void get_opts(int argc, char **argv, struct options_t *opts);
void get_opts(int argc, char **argv, options_t *opts)
{
  if (argc == 1)
  {
    std::cout << "Usage:\n"
              << "\t--in or -i <file_path>\n"
              << "\t--out or -o <file_path>\n"
              << "\t--steps or -s <steps>\n"
              << "\t--theta or -t <theta>\n"
              << "\t--timestep or -d <timestep>\n"
              << "\t--visualization or -v <visualization>\n";
    std::exit(0);
  }

  opts->in_file = nullptr;
  opts->out_file = nullptr;
  opts->steps = 100;
  opts->theta = 0.7;
  opts->timestep = 0.01;
  opts->visualization = false;

  static struct option l_opts[] = {
      {"in", required_argument, nullptr, 'i'},
      {"out", required_argument, nullptr, 'o'},
      {"steps", required_argument, nullptr, 's'},
      {"theta", required_argument, nullptr, 't'},
      {"timestep", required_argument, nullptr, 'd'},
      {"visualization", required_argument, nullptr, 'v'},
      {0, 0, 0, 0} // terminator
  };

  int opt, ind;

  while ((opt = getopt_long(argc, argv, "i:o:s:t:d:v:", l_opts, &ind)) != -1)
  {
    switch (opt)
    {
    case 'i':
      opts->in_file = optarg;
      break;

    case 'o':
      opts->out_file = optarg;
      break;

    case 's':
      opts->steps = std::atoi(optarg);
      break;

    case 't':
      opts->theta = std::atof(optarg);
      break;

    case 'd':
      opts->timestep = std::atof(optarg);
      break;

    case 'v':
      opts->visualization = optarg;
      break;

    case '?':

      std::cerr << "unkown arg\n";
      std::exit(1);

    default:
      break;
    }
  }
}

struct Particle
{
  int index;
  double x;
  double y;
  double mass;
  double vx;
  double vy;
};

struct Node
{
  double x_min, x_max;
  double y_min, y_max;

  Node *upper_left;
  Node *upper_right;
  Node *lower_left;
  Node *lower_right;

  double total_mass;
  double center_x;
  double center_y;

  Particle *particle;

  Node(double xmin, double xmax, double ymin, double ymax)
      : x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax),
        upper_left(nullptr), upper_right(nullptr), lower_left(nullptr), lower_right(nullptr),
        particle(nullptr),
        total_mass(0), center_x(0), center_y(0)
  {
  }
};
void write_file(const char *filename, const std::vector<Particle> &particles)
{
  std::ofstream outfile(filename);
  if (!outfile)
  {
    std::cerr << "Error: could not open output file " << filename << "\n";
    return;
  }

  // Optionally write number of particles
  outfile << particles.size() << "\n";

  for (const auto &p : particles)
  {
    outfile << p.index << " "
            << p.x << " "
            << p.y << " "
            << p.mass << " "
            << p.vx << " "
            << p.vy << "\n";
  }

  outfile.close();
}

int read_file(struct options_t *args,
              std::vector<Particle> &particles)
{
  std::ifstream infile(args->in_file);
  if (!infile)
  {
    std::cerr << "Error opening file\n";
    return 1;
  }

  int n;
  infile >> n;

  particles.resize(n);

  for (int i = 0; i < n; i++)
  {
    Particle b;
    infile >> b.index >> b.x >> b.y >> b.mass >> b.vx >> b.vy;

    if (!infile)
    {
      std::cerr << "Error reading body " << i << std::endl;
      return 1;
    }

    particles[i] = b;
  }

  return 0;
}

double distance(Particle a, Particle b)
{
  double x_dist = a.x - b.x;
  double y_dist = a.y - b.y;
  double dist = sqrt(pow(x_dist, 2) + pow(y_dist, 2));
  return std::max(dist, rlimit);
}

double distance(Particle *a, Particle *b)
{
  double x_dist = a->x - b->x;
  double y_dist = a->y - b->y;
  double dist = sqrt(pow(x_dist, 2) + pow(y_dist, 2));
  return std::max(dist, rlimit);
}

double distanceP2N(Particle *a, Node *b)
{
  double x_dist = a->x - b->center_x;
  double y_dist = a->y - b->center_y;
  double dist = sqrt(pow(x_dist, 2) + pow(y_dist, 2));
  return std::max(dist, rlimit);
}

double calcForceX(Particle a, Particle b)
{
  double dist = distance(a, b);
  double dist_x = a.x - b.x;
  return -(G * a.mass * b.mass * dist_x) / pow(dist, 3);
}

double calcForceY(Particle a, Particle b)
{
  double dist = distance(a, b);
  double dist_y = a.y - b.y;
  return -(G * a.mass * b.mass * dist_y) / pow(dist, 3);
}

double calcForceX(Particle *a, Node *b)
{
  double dist = distanceP2N(a, b);
  double dist_x = a->x - b->center_x;
  return -(G * a->mass * b->total_mass * dist_x) / pow(dist, 3);
}

double calcForceY(Particle *a, Node *b)
{
  double dist = distanceP2N(a, b);
  double dist_y = a->y - b->center_y;
  return -(G * a->mass * b->total_mass * dist_y) / pow(dist, 3);
}

double calcForceX(Particle *a, Particle *b)
{
  double dist = distance(a, b);
  double dist_x = a->x - b->x;
  return -(G * a->mass * b->mass * dist_x) / pow(dist, 3);
}

double calcForceY(Particle *a, Particle *b)
{
  double dist = distance(a, b);
  double dist_y = a->y - b->y;
  return -(G * a->mass * b->mass * dist_y) / pow(dist, 3);
}

// void drawParticle2D(double x_window, double y_window,
//                     double radius,
//                     float *colors)
// {
//   int k = 0;
//   float angle = 0.0f;
//   glBegin(GL_TRIANGLE_FAN);
//   glColor3f(colors[0], colors[1], colors[2]);
//   glVertex2f(x_window, y_window);
//   for (k = 0; k < 20; k++)
//   {
//     angle = (float)(k) / 19 * 2 * 3.141592;
//     glVertex2f(x_window + radius * cos(angle),
//                y_window + radius * sin(angle));
//   }
//   glEnd();
// }

void insertNode(Node *node, Particle *particle)
{
  double vert_line = (node->x_max + node->x_min) / 2.0;
  double horz_line = (node->y_max + node->y_min) / 2.0;

  bool left = particle->x < vert_line;
  bool down = particle->y < horz_line;

  // if there was a particle and we are trying to insert, we must subdivide and reinsert
  if (node->particle != nullptr)
  {
    Particle *existing = node->particle;
    node->particle = nullptr;

    insertNode(node, existing);
  }

  if (left)
  {
    // lower left
    if (down)
    {
      if (node->lower_left)
      {
        insertNode(node->lower_left, particle);
      }
      else
      {
        node->lower_left = new Node(node->x_min, vert_line, node->y_min, horz_line);
        node->lower_left->particle = particle;
      }
    }
    // upper left
    else
    {
      if (node->upper_left)
      {
        insertNode(node->upper_left, particle);
      }
      else
      {
        node->upper_left = new Node(node->x_min, vert_line, horz_line, node->y_max);
        node->upper_left->particle = particle;
      }
    }
  }
  else
  {
    // lower right
    if (down)
    {
      if (node->lower_right)
      {
        insertNode(node->lower_right, particle);
      }
      else
      {
        node->lower_right = new Node(vert_line, node->x_max, node->y_min, horz_line);
        node->lower_right->particle = particle;
      }
    }
    // upper right
    else
    {
      if (node->upper_right)
      {
        insertNode(node->upper_right, particle);
      }
      else
      {
        node->upper_right = new Node(vert_line, node->x_max, horz_line, node->y_max);
        node->upper_right->particle = particle;
      }
    }
  }
}

double centerOfMass(double pos_a, double pos_b, double mass_a, double mass_b)
{
  return ((pos_a * mass_a) + (pos_b * mass_b)) / (mass_a + mass_b);
}

void calculateCenters(Node *node)
{

  if (node->lower_left)
  {
    if (node->lower_left->particle)
    {
      node->center_x = centerOfMass(node->lower_left->particle->x, node->center_x, node->lower_left->particle->mass, node->total_mass);
      node->center_y = centerOfMass(node->lower_left->particle->y, node->center_y, node->lower_left->particle->mass, node->total_mass);
      node->total_mass += node->lower_left->particle->mass;
    }
    else
    {
      calculateCenters(node->lower_left);
      node->center_x = centerOfMass(node->lower_left->center_x, node->center_x, node->lower_left->total_mass, node->total_mass);
      node->center_y = centerOfMass(node->lower_left->center_y, node->center_y, node->lower_left->total_mass, node->total_mass);
      node->total_mass += node->lower_left->total_mass;
    }
  }
  if (node->lower_right)
  {

    if (node->lower_right->particle)
    {
      node->center_x = centerOfMass(node->lower_right->particle->x, node->center_x, node->lower_right->particle->mass, node->total_mass);
      node->center_y = centerOfMass(node->lower_right->particle->y, node->center_y, node->lower_right->particle->mass, node->total_mass);
      node->total_mass += node->lower_right->particle->mass;
    }
    else
    {
      calculateCenters(node->lower_right);
      node->center_x = centerOfMass(node->lower_right->center_x, node->center_x, node->lower_right->total_mass, node->total_mass);
      node->center_y = centerOfMass(node->lower_right->center_y, node->center_y, node->lower_right->total_mass, node->total_mass);
      node->total_mass += node->lower_right->total_mass;
    }
  }
  if (node->upper_right)
  {

    if (node->upper_right->particle)
    {
      node->center_x = centerOfMass(node->upper_right->particle->x, node->center_x, node->upper_right->particle->mass, node->total_mass);
      node->center_y = centerOfMass(node->upper_right->particle->y, node->center_y, node->upper_right->particle->mass, node->total_mass);
      node->total_mass += node->upper_right->particle->mass;
    }
    else
    {
      calculateCenters(node->upper_right);
      node->center_x = centerOfMass(node->upper_right->center_x, node->center_x, node->upper_right->total_mass, node->total_mass);
      node->center_y = centerOfMass(node->upper_right->center_y, node->center_y, node->upper_right->total_mass, node->total_mass);
      node->total_mass += node->upper_right->total_mass;
    }
  }

  if (node->upper_left)
  {
    if (node->upper_left->particle)
    {
      node->center_x = centerOfMass(node->upper_left->particle->x, node->center_x, node->upper_left->particle->mass, node->total_mass);
      node->center_y = centerOfMass(node->upper_left->particle->y, node->center_y, node->upper_left->particle->mass, node->total_mass);
      node->total_mass += node->upper_left->particle->mass;
    }
    else
    {
      calculateCenters(node->upper_left);
      node->center_x = centerOfMass(node->upper_left->center_x, node->center_x, node->upper_left->total_mass, node->total_mass);
      node->center_y = centerOfMass(node->upper_left->center_y, node->center_y, node->upper_left->total_mass, node->total_mass);
      node->total_mass += node->upper_left->total_mass;
    }
  }
}
void computeNaive(std::vector<Particle> &particles, double dt)
{
  std::vector<Particle> next_particles = particles;

  for (int i = 0; i < (int)particles.size(); i++)
  {
    const Particle &current_particle = particles[i];

    // If this particle is already dead, skip
    if (current_particle.mass < 0)
      continue;

    if (current_particle.x < 0 || current_particle.x > MAX_X ||
        current_particle.y < 0 || current_particle.y > MAX_Y)
    {
      next_particles[i].mass = -1;
      continue;
    }

    double total_x_f = 0.0;
    double total_y_f = 0.0;

    for (int j = 0; j < (int)particles.size(); j++)
    {
      const Particle &other = particles[j];

      if (other.mass < 0)
        continue;

      if (other.x < 0 || other.x > MAX_X ||
          other.y < 0 || other.y > MAX_Y)
      {
        next_particles[j].mass = -1;
        continue;
      }

      if (i == j)
        continue;

      double force_x = calcForceX(current_particle, other);
      double force_y = calcForceY(current_particle, other);
      total_x_f += force_x;
      total_y_f += force_y;
    }

    double acc_x = total_x_f / current_particle.mass;
    double acc_y = total_y_f / current_particle.mass;

    double dt2 = dt * dt;
    double new_pos_x = current_particle.x + current_particle.vx * dt + 0.5 * acc_x * dt2;
    double new_pos_y = current_particle.y + current_particle.vy * dt + 0.5 * acc_y * dt2;

    double new_vx = current_particle.vx + acc_x * dt;
    double new_vy = current_particle.vy + acc_y * dt;

    next_particles[i].x = new_pos_x;
    next_particles[i].y = new_pos_y;
    next_particles[i].vx = new_vx;
    next_particles[i].vy = new_vy;
  }

  particles = next_particles;
}

static inline void addForce(
    double dx,
    double dy,
    double m1,
    double m2,
    double *fx,
    double *fy)
{
  double r2 = dx * dx + dy * dy;
  double r = sqrt(r2);
  double dist = std::max(rlimit, r);
  double inv_r3 = 1.0 / (dist * dist * dist);

  double F = G * m1 * m2 * inv_r3;

  *fx += F * dx;
  *fy += F * dy;
}

MPI_Datatype MPI_PARTICLE;

struct ParticleUpdate
{
  double x;
  double y;
  double vx;
  double vy;
};

void create_particle_update_type()
{
  const int nitems = 4;

  int blocklengths[4] = {1, 1, 1, 1};
  MPI_Datatype types[4] = {
      MPI_DOUBLE, // x
      MPI_DOUBLE, // y
      MPI_DOUBLE, // vx
      MPI_DOUBLE, // vy

  };

  MPI_Aint offsets[4];
  offsets[0] = offsetof(ParticleUpdate, x);
  offsets[1] = offsetof(ParticleUpdate, y);
  offsets[2] = offsetof(ParticleUpdate, vx);
  offsets[3] = offsetof(ParticleUpdate, vy);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);
}

void calculateNodeForce(Particle *p, Node *n, double *total_x_f, double *total_y_f, double theta2)
{
  if (n->total_mass < 0.0)
    return;

  if (n->particle)
  {
    if (n->particle->index == p->index)
      return;
    if (n->particle->mass < 0)
      return;

    double dx = n->particle->x - p->x;
    double dy = n->particle->y - p->y;
    addForce(dx, dy, p->mass, n->particle->mass, total_x_f, total_y_f);
    return;
  }

  double dx = n->center_x - p->x;
  double dy = n->center_y - p->y;
  double r2 = dx * dx + dy * dy;
  double l = n->x_max - n->x_min;

  if (r2 > 0.0 && (l * l) < theta2 * r2)
  {
    addForce(dx, dy, p->mass, n->total_mass, total_x_f, total_y_f);
    return;
  }

  if (n->lower_left)
    calculateNodeForce(p, n->lower_left, total_x_f, total_y_f, theta2);
  if (n->lower_right)
    calculateNodeForce(p, n->lower_right, total_x_f, total_y_f, theta2);
  if (n->upper_left)
    calculateNodeForce(p, n->upper_left, total_x_f, total_y_f, theta2);
  if (n->upper_right)
    calculateNodeForce(p, n->upper_right, total_x_f, total_y_f, theta2);
}
void deleteTree(Node *node)
{
  if (!node)
    return;
  deleteTree(node->lower_left);
  deleteTree(node->lower_right);
  deleteTree(node->upper_left);
  deleteTree(node->upper_right);
  delete node;
}
void computeBarnesHunt(std::vector<Particle> &particles, double dt, double theta)
{
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // build the tree
  Node *root = new Node(0, MAX_X, 0, MAX_Y);

  // seperate nodes into quadrants
  for (int i = 0; i < particles.size(); i++)
  {
    if (particles[i].x < 0 || particles[i].x > MAX_X ||
        particles[i].y < 0 || particles[i].y > MAX_Y)
    {
      particles[i].mass = -1;
      continue;
    }
    insertNode(root, &particles[i]);
  }

  // find center of mass for every tree node
  calculateCenters(root);

  std::vector<Particle> next_particles = particles;

  // calculate for force for every particle based on tree nodes
  const std::size_t n = particles.size();
  const double dt2 = dt * dt;
  const double half_dt2 = 0.5 * dt2;
  const double theta2 = theta * theta;

  int slice = (n + world_size - 1) / world_size;
  int start = slice * world_rank;
  int end = std::min(slice * (world_rank + 1), (int)n);

  for (int i = start; i < end; i++)
  {
    if (particles[i].mass < 0 ||
        particles[i].x < 0 || particles[i].x > MAX_X ||
        particles[i].y < 0 || particles[i].y > MAX_Y)
    {
      continue;
    }

    const Particle &current_particle = particles[i];

    double force_x = 0;
    double force_y = 0;

    calculateNodeForce(&particles[i], root, &force_x, &force_y, theta2);

    double acc_x = force_x / current_particle.mass;
    double acc_y = force_y / current_particle.mass;

    double new_pos_x = current_particle.x + (current_particle.vx * dt) + (acc_x * half_dt2);
    double new_pos_y = current_particle.y + (current_particle.vy * dt) + (acc_y * half_dt2);

    double new_vx = current_particle.vx + acc_x * dt;
    double new_vy = current_particle.vy + acc_y * dt;

    next_particles[i].x = new_pos_x;
    next_particles[i].y = new_pos_y;
    next_particles[i].vx = new_vx;
    next_particles[i].vy = new_vy;
  }

  if (world_size > 1)
  {
    if (world_rank == 0)
    {
      for (int i = 1; i < world_size; i++)
      {
        std::vector<ParticleUpdate> particle_slice(slice);
        MPI_Recv(particle_slice.data(), slice, MPI_PARTICLE, i, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int j = 0; j < slice; j++)
        {
          int idx = j + i * slice;
          if (idx >= n)
            break;

          next_particles[idx].x = particle_slice[j].x;
          next_particles[idx].y = particle_slice[j].y;
          next_particles[idx].vx = particle_slice[j].vx;
          next_particles[idx].vy = particle_slice[j].vy;
        }
      }
    }
    else
    {
      std::vector<ParticleUpdate> particle_slice(slice);

      for (int i = 0; i < slice; i++)
      {
        int idx = i + world_rank * slice;
        if (idx >= n)
          break;

        particle_slice[i].x = next_particles[idx].x;
        particle_slice[i].y = next_particles[idx].y;
        particle_slice[i].vx = next_particles[idx].vx;
        particle_slice[i].vy = next_particles[idx].vy;
      }

      MPI_Send(particle_slice.data(), slice, MPI_PARTICLE, 0, 0, MPI_COMM_WORLD);
    }

    if (world_rank == 0)
    {
      std::vector<ParticleUpdate> all_updates(n);
      for (int j = 0; j < n; j++)
      {
        all_updates[j].x = next_particles[j].x;
        all_updates[j].y = next_particles[j].y;
        all_updates[j].vx = next_particles[j].vx;
        all_updates[j].vy = next_particles[j].vy;
      }

      for (int i = 1; i < world_size; i++)
      {
        MPI_Send(all_updates.data(), n, MPI_PARTICLE, i, 0, MPI_COMM_WORLD);
      }
    }
    else
    {
      std::vector<ParticleUpdate> all_updates(n);
      MPI_Recv(all_updates.data(), n, MPI_PARTICLE, 0, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for (int j = 0; j < n; j++)
      {
        next_particles[j].x = all_updates[j].x;
        next_particles[j].y = all_updates[j].y;
        next_particles[j].vx = all_updates[j].vx;
        next_particles[j].vy = all_updates[j].vy;
      }
    }
  }

  particles = next_particles;
  deleteTree(root);
}

int main(int argc, char **argv)
{
  struct options_t opts;
  get_opts(argc, argv, &opts);
  // GLFWwindow *window;
  MPI_Init(NULL, NULL);

  create_particle_update_type();
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // if (opts.visualization && world_rank == 0)
  // {
  //   /* OpenGL window dims */
  //   int width = 600, height = 600;
  //   if (!glfwInit())
  //   {
  //     fprintf(stderr, "Failed to initialize GLFW\n");
  //     return -1;
  //   }
  //   // Open a window and create its OpenGL context
  //   window = glfwCreateWindow(width, height, "Simulation", NULL, NULL);
  //   if (window == NULL)
  //   {
  //     fprintf(stderr, "Failed to open GLFW window.\n");
  //     glfwTerminate();
  //     return -1;
  //   }
  //   glfwMakeContextCurrent(window); // Initialize GLEW
  //   if (glewInit() != GLEW_OK)
  //   {
  //     fprintf(stderr, "Failed to initialize GLEW\n");
  //     return -1;
  //   }
  //   // Ensure we can capture the escape key being pressed below
  //   glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
  // }

  std::vector<Particle> particles;
  if (read_file(&opts, particles) != 0)
  {
    std::cerr << "Failed to read file\n";
  }

  double start = MPI_Wtime();

  for (int t = 0; t < opts.steps; t++)
  {
    // if (opts.visualization && world_rank == 0)
    // {
    //   glClear(GL_COLOR_BUFFER_BIT);

    //   for (int p = 0; p < particles.size(); p++)
    //   {
    //     if (particles[p].mass < 0)
    //       continue;

    //     double x_window = 2 * particles[p].x / MAX_X - 1;
    //     double y_window = 2 * particles[p].y / MAX_Y - 1;

    //     float colors[3] = {0.01f, 0.5f, 0.99f};
    //     drawParticle2D(x_window, y_window, 0.01, colors);
    //   }

    //   glfwSwapBuffers(window);
    //   glfwPollEvents();
    // }

    computeBarnesHunt(particles, opts.timestep, opts.theta);
    // computeNaive(particles);
  }
  double end = MPI_Wtime();
  double elapsed = end - start;

  if (world_rank == 0)
  {
    std::cout << elapsed << "\n";
  }

  if (opts.out_file && world_rank == 0)
  {
    write_file(opts.out_file, particles);
  }

  MPI_Type_free(&MPI_PARTICLE);

  MPI_Finalize();
}