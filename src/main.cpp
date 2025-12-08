
#include <mpi.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <GL/glu.h>
#include <GL/glut.h>

const double G = 0.0001;
const double rlimit = 0.03;
const double dt = 0.01;
const double theta = 0.05;

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

  // Set some sane defaults (optional)
  opts->in_file = nullptr;
  opts->out_file = nullptr;
  opts->steps = 0;
  opts->theta = 0.0;
  opts->timestep = 0.0;
  opts->visualization = false;

  // Long options table
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
  // Short options: each option that expects an argument gets a ':' after it
  // i:o:s:t:d:v:
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
      // Unknown option or missing argument.
      // getopt_long already printed an error message.
      std::cerr << "Use --help or run with no args to see usage.\n";
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

  // Resize vector to hold n particles
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

void drawParticle2D(double x_window, double y_window,
                    double radius,
                    float *colors)
{
  int k = 0;
  float angle = 0.0f;
  glBegin(GL_TRIANGLE_FAN);
  glColor3f(colors[0], colors[1], colors[2]);
  glVertex2f(x_window, y_window);
  for (k = 0; k < 20; k++)
  {
    angle = (float)(k) / 19 * 2 * 3.141592;
    glVertex2f(x_window + radius * cos(angle),
               y_window + radius * sin(angle));
  }
  glEnd();
}

// void drawOctreeBounds2D(Particle node)
// {
//   int i;
//   if (!node || node is external)
//     return;
//   glBegin(GL_LINES);
//   // set the color of lines to be white
//   glColor3f(1.0f, 1.0f, 1.0f);
//   // specify the start point's coordinates
//   glVertex2f(h_start_x, h_start_y);
//   // specify the end point's coordinates
//   glVertex2f(h_end_x, h_end_y);
//   // do the same for verticle line
//   glVertex2f(v_start_x, h_end_y);
//   glVertex2f(v_end_x, v_end_y);
//   glEnd();
//     for each
//       subtree of node
//           drawOctreeBounds2D(subtree);
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
    node->lower_left = new Node{node->x_min, vert_line, node->y_min, horz_line};
    node->lower_right = new Node{vert_line, node->x_max, node->y_min, horz_line};
    node->upper_left = new Node{node->x_min, vert_line, horz_line, node->y_max};
    node->upper_right = new Node{vert_line, node->x_max, horz_line, node->y_max};
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

void computeNaive(std::vector<Particle> &particles)
{

  std::vector<Particle> next_particles = particles;

  for (int i = 0; i < particles.size(); i++)
  {
    const Particle &current_particle = particles[i];
    double total_x_f = 0;
    double total_y_f = 0;

    for (int j = 0; j < particles.size(); j++)
    {
      if (i == j)
        continue;

      // calculate the total force from every partical, then figure out acceleration

      double force_x = calcForceX(current_particle, particles[j]);
      double force_y = calcForceY(current_particle, particles[j]);
      total_x_f += force_x;
      total_y_f += force_y;
    }

    double acc_x = total_x_f / current_particle.mass;
    double acc_y = total_y_f / current_particle.mass;

    double new_pos_x = current_particle.x + (current_particle.vx * dt) + (0.5 * acc_x * pow(dt, 2));
    double new_pos_y = current_particle.y + (current_particle.vy * dt) + (0.5 * acc_y * pow(dt, 2));

    double new_vx = current_particle.vx + acc_x * dt;
    double new_vy = current_particle.vy + acc_y * dt;

    next_particles[i].x = new_pos_x;
    next_particles[i].y = new_pos_y;
    next_particles[i].vx = new_vx;
    next_particles[i].vy = new_vy;
  }

  particles = next_particles;
}

void calculateNodeForce(Particle *p, Node *n, double *total_x_f, double *total_y_f)
{

  // if we got to a particle, lets just add
  if (n->particle)
  {

    double force_x = calcForceX(p, n->particle);
    double force_y = calcForceY(p, n->particle);

    *total_x_f += force_x;
    *total_y_f += force_y;

    return;
  }

  double l = n->x_max - n->x_min;
  double d = distanceP2N(p, n);

  if ((l / d) < theta)
  {

    double force_x = calcForceX(p, n);
    double force_y = calcForceY(p, n);
    *total_x_f += force_x;
    *total_y_f += force_y;
  }
  else
  {
    if (n->lower_left)
    {
      calculateNodeForce(p, n->lower_left, total_x_f, total_y_f);
    }
    if (n->lower_right)
    {
      calculateNodeForce(p, n->lower_right, total_x_f, total_y_f);
    }
    if (n->upper_left)
    {
      calculateNodeForce(p, n->upper_left, total_x_f, total_y_f);
    }
    if (n->upper_right)
    {
      calculateNodeForce(p, n->upper_right, total_x_f, total_y_f);
    }
  }
}

void computeBarnesHunt(std::vector<Particle> &particles)
{

  // build the tree
  Node *root = new Node(0, MAX_X, 0, MAX_Y);

  double horz_line = MAX_X / 2.0;
  double vert_line = MAX_Y / 2.0;

  // seperate nodes into quadrants
  for (int i = 0; i < particles.size(); i++)
  {
    if (particles[i].x < 0 || particles[i].x > MAX_X || particles[i].y < 0 || particles[i].y > MAX_Y)
    {
      particles[i].mass = 0;
      continue;
    }
    insertNode(root, &particles[i]);
  }

  // find center of mass for every tree node
  // sum masses and add up x and y positions
  calculateCenters(root);

  std::vector<Particle> next_particles = particles;

  // caclulate for force for every particle based on tree nodes
  for (int i = 0; i < particles.size(); i++)
  {
    if (particles[i].x < 0 || particles[i].x > MAX_X || particles[i].y < 0 || particles[i].y > MAX_Y)
    {
      continue;
    }
    const Particle &current_particle = particles[i];

    double force_x = 0;
    double force_y = 0;

    calculateNodeForce(&particles[i], root, &force_x, &force_y);

    double acc_x = force_x / current_particle.mass;
    double acc_y = force_y / current_particle.mass;

    double new_pos_x = current_particle.x + (current_particle.vx * dt) + (0.5 * acc_x * pow(dt, 2));
    double new_pos_y = current_particle.y + (current_particle.vy * dt) + (0.5 * acc_y * pow(dt, 2));

    double new_vx = current_particle.vx + acc_x * dt;
    double new_vy = current_particle.vy + acc_y * dt;

    next_particles[i].x = new_pos_x;
    next_particles[i].y = new_pos_y;
    next_particles[i].vx = new_vx;
    next_particles[i].vy = new_vy;
  }
  particles = next_particles;
}

int main(int argc, char **argv)
{

  /* OpenGL window dims */
  int width = 600, height = 600;
  GLFWwindow *window;
  if (!glfwInit())
  {
    fprintf(stderr, "Failed to initialize GLFW\n");
    return -1;
  }
  // Open a window and create its OpenGL context
  window = glfwCreateWindow(width, height, "Simulation", NULL, NULL);
  if (window == NULL)
  {
    fprintf(stderr, "Failed to open GLFW window.\n");
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(window); // Initialize GLEW
  if (glewInit() != GLEW_OK)
  {
    fprintf(stderr, "Failed to initialize GLEW\n");
    return -1;
  }
  // Ensure we can capture the escape key being pressed below
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

  struct options_t opts;
  get_opts(argc, argv, &opts);

  std::vector<Particle> particles;
  if (read_file(&opts, particles) != 0)
  {
    std::cerr << "Failed to read file\n";
  }
  for (int i = 0; i < particles.size(); i++)
  {
    std::cout << "p" << i << ": x=" << particles[i].x
              << ", y=" << particles[i].y << "\n";
  }

  // Initialize the MPI environment. The two arguments to MPI Init are not
  // currently used by MPI implementations, but are there in case future
  // implementations might need the arguments.
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d out of %d processors\n",
         processor_name, world_rank, world_size);

  // Finalize the MPI environment. No more MPI calls can be made after this

  // testing most basic n * n implementation

  // lets do for 100 time steps
  for (int t = 0; t < 1000000; t++)
  {

    glClear(GL_COLOR_BUFFER_BIT);

    // drawOctreeBounds2D(node);
    for (int p = 0; p < particles.size(); p++)
    {
      if (particles[p].mass == 0)
        continue;
        
      double x_window = 2 * particles[p].x / MAX_X - 1;
      double y_window = 2 * particles[p].y / MAX_Y - 1;

      float colors[3] = {0.01f, 0.5f, 0.99f};
      drawParticle2D(x_window, y_window, 0.03, colors);
    }

    // Swap buffers
    glfwSwapBuffers(window);
    glfwPollEvents();

    computeBarnesHunt(particles);
  }

  std::cout << "final post\n";

  for (int i = 0; i < particles.size(); i++)
  {
    std::cout << "p" << i << ": x=" << particles[i].x
              << ", y=" << particles[i].y << "\n";
  }

  MPI_Finalize();
}