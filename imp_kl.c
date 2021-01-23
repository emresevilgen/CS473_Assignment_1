#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <string.h>
#include "mmio.h"

// Adjacency List
struct node
{
    int vertex;
    struct node *next;
};

struct Graph
{
    int numVertices;
    struct node **adjLists;
};

struct node *createNode(int v)
{
    struct node *newNode = malloc(sizeof(struct node));
    newNode->vertex = v;
    newNode->next = NULL;
    return newNode;
}

struct Graph *createAGraph(int vertices)
{
    struct Graph *graph = malloc(sizeof(struct Graph));
    graph->numVertices = vertices;

    graph->adjLists = malloc(vertices * sizeof(struct node *));

    int i;
    for (i = 0; i < vertices; i++)
        graph->adjLists[i] = NULL;

    return graph;
}

void addEdge(struct Graph *graph, int s, int d)
{
    struct node *newNode = createNode(d);
    newNode->next = graph->adjLists[s];
    graph->adjLists[s] = newNode;

    newNode = createNode(s);
    newNode->next = graph->adjLists[d];
    graph->adjLists[d] = newNode;
}

// Heap
struct heap_node
{
    int gain;
    int index;
};

void swap_indexes(int *a, int *b)
{
    int temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

void swap(struct heap_node *a, struct heap_node *b, int *indexes)
{
    swap_indexes(&indexes[(*a).index], &indexes[(*b).index]);

    struct heap_node temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

int get_parent(struct heap_node *A, int size, int index)
{
    if ((index > 0) && (index < size))
    {
        return (index - 1) / 2;
    }
    return -1;
}

int get_left_child(struct heap_node *A, int size, int index)
{
    if (((2 * index + 1) < size) && (index >= 0))
        return 2 * index + 1;
    return -1;
}

int get_right_child(struct heap_node *A, int size, int index)
{
    if ((((2 * index) + 2) < size) && (index >= 0))
        return (2 * index) + 2;
    return -1;
}

void heapify(struct heap_node *A, int *A_indexes, int size, int index)
{
    int left_child_index = get_left_child(A, size, index);
    int right_child_index = get_right_child(A, size, index);

    // Finding largest of current, left child and right child
    int largest = index;

    if ((left_child_index < size) && (left_child_index > 0))
    {
        if (A[left_child_index].gain > A[largest].gain)
        {
            largest = left_child_index;
        }
    }

    if ((right_child_index < size && (right_child_index > 0)))
    {
        if (A[right_child_index].gain > A[largest].gain)
        {
            largest = right_child_index;
        }
    }

    // Largest is not the current
    if (largest != index)
    {
        swap(&A[index], &A[largest], A_indexes);
        heapify(A, A_indexes, size, largest);
    }
}

void build_heap(struct heap_node *A, int *A_indexes, int size)
{
    int i;
    for (i = size; i >= 0; i--)
    {
        heapify(A, A_indexes, size, i);
    }
}

struct heap_node extract_max(struct heap_node *A, int *A_indexes, int *sizePtr)
{
    struct heap_node root = A[0];
    *sizePtr = *sizePtr - 1;
    swap(&A[0], &A[*sizePtr], A_indexes);
    heapify(A, A_indexes, *sizePtr, 0);
    return root;
}

void increase_key(struct heap_node *A, int *A_indexes, int size, int index, int key)
{
    A[index].gain = key;
    while ((index > 0) && (A[get_parent(A, size, index)].gain < A[index].gain))
    {
        swap(&A[index], &A[get_parent(A, size, index)], A_indexes);
        index = get_parent(A, size, index);
    }
}

void decrease_key(struct heap_node *A, int *A_indexes, int size, int index, int key)
{
    A[index].gain = key;
    heapify(A, A_indexes, size, index);
}

void print_heap(struct heap_node *A, int size)
{
    int i;
    int j = 0;
    for (i = 0; i < size; i++)
    {
        if (i == (1 << j) - 1)
        {
            printf("\n");
            j++;
        }
        printf("%d ", A[i].gain);
    }
    printf("\n");
}

struct location
{
    int a_or_b;
    int index;
};

// Calculate cut value
int calculate_cut_value(struct Graph *graph, int *A, int size_of_A, int *B, int size_of_B, struct location *locations)
{
    int i;
    int cut_value = 0;
    struct node *cur;
    struct node *next;

    for (i = 0; i < size_of_A; i++)
    {
        cur = (graph->adjLists)[A[i]];
        while (cur != NULL)
        {
            if (locations[cur->vertex].a_or_b == 1)
            {
                cut_value++;
            }
            next = cur->next;
            cur = next;
        }
    }
    return cut_value;
}

// Read file
int readMatrixFile(FILE *file, struct Graph **graphPtr, int *no_of_verticesPtr)
{
    MM_typecode matcode;

    if (mm_read_banner(file, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        return 1;
    }

    int ret_code;
    int M, N, nz;

    if ((ret_code = mm_read_mtx_crd_size(file, &M, &N, &nz)) != 0)
        return 1;

    *graphPtr = createAGraph(M);
    struct Graph *graph = *graphPtr;

    int i;
    int src;
    int dest;

    for (i = 0; i < nz; i++)
    {
        int ret = fscanf(file, "%d %d\n", &src, &dest);
        addEdge(graph, src - 1, dest - 1);
    }

    *no_of_verticesPtr = M;

    if (file != stdin)
        fclose(file);

    return 0;
}

// Simpler KL implementations
int simpler_KL(struct Graph *graph, int *initial_cut_sizePtr, int *final_cut_sizePtr, double *time_spentPtr)
{

    // Split V into two balanced disjoint sets A and B
    int no_of_vertices = graph->numVertices;
    struct location *locations = malloc(no_of_vertices * sizeof(struct location));

    int size_of_A = no_of_vertices / 2;
    int *A = malloc(size_of_A * sizeof(int));
    int i;
    int j;

    for (i = 0; i < size_of_A; i++)
    {
        A[i] = i;
        locations[i].a_or_b = 0;
        locations[i].index = i;
    }

    int size_of_B = no_of_vertices - no_of_vertices / 2;
    int *B = malloc(size_of_B * sizeof(int));

    for (j = 0; j < size_of_B; j++)
    {
        B[j] = i;
        locations[i].a_or_b = 1;
        locations[i].index = j;
        i++;
    }

    // Calculate initial cut value
    int cut_value = calculate_cut_value(graph, A, size_of_A, B, size_of_B, locations);
    *initial_cut_sizePtr = cut_value;

    int gmax = 0;

    // Create lock array
    int *isLocked = malloc(no_of_vertices * sizeof(int));

    clock_t begin = clock();

    struct node *cur;
    struct node *next;

    // Loop mallocs
    int *gain_A = malloc(size_of_A * sizeof(int));
    int *gain_B = malloc(size_of_B * sizeof(int));

    int *Lg = malloc(size_of_A * sizeof(int));
    int *La = malloc(size_of_A * sizeof(int));
    int *Lb = malloc(size_of_A * sizeof(int));

    do
    {
        // Initialize lock array
        for (i = 0; i < no_of_vertices; i++)
            isLocked[i] = 0;

        // Compute D values (move gain values) for all A and B

        for (i = 0; i < size_of_A; i++)
            gain_A[i] = 0;
        for (i = 0; i < size_of_B; i++)
            gain_B[i] = 0;

        for (i = 0; i < size_of_A; i++)
        {
            cur = (graph->adjLists)[A[i]];
            while (cur != NULL)
            {
                struct location loc = locations[cur->vertex];

                if (loc.a_or_b == 0)
                {
                    gain_A[i] = gain_A[i] - 1;
                }
                else
                {
                    gain_A[i] = gain_A[i] + 1;
                }
                next = cur->next;
                cur = next;
            }
        }

        for (i = 0; i < size_of_B; i++)
        {
            cur = (graph->adjLists)[B[i]];
            while (cur != NULL)
            {
                struct location loc = locations[cur->vertex];

                if (loc.a_or_b == 1)
                {
                    gain_B[i] = gain_B[i] - 1;
                }
                else
                {
                    gain_B[i] = gain_B[i] + 1;
                }
                next = cur->next;
                cur = next;
            }
        }

        // Let Lg, La, Lb be empty lists and i = 0

        int i_L = 0;

        int n;
        for (n = 0; n < size_of_A; n++)
        {
            // Find a and b such that the exchange gain g = Da + Db is maximized.
            int max_A_index = -1;
            int da = INT_MIN;
            for (j = 0; j < size_of_A; j++)
            {
                if (isLocked[A[j]] == 0 && da < gain_A[j])
                {
                    max_A_index = j;
                    da = gain_A[j];
                }
            }

            if (max_A_index == -1)
                break;

            int max_B_index = -1;
            int db = INT_MIN;
            for (j = 0; j < size_of_B; j++)
            {
                if (isLocked[B[j]] == 0 && db < gain_B[j])
                {
                    max_B_index = j;
                    db = gain_B[j];
                }
            }

            if (max_B_index == -1)
                break;

            int gain = gain_A[max_A_index] + gain_B[max_B_index];
            int a_index_A = max_A_index;
            int b_index_B = max_B_index;

            max_A_index = A[max_A_index];
            max_B_index = B[max_B_index];

            cur = (graph->adjLists)[max_A_index];
            while (cur != NULL)
            {
                if (cur->vertex == max_B_index)
                {
                    gain = gain - 2;
                    break;
                }
                next = cur->next;
                cur = next;
            }

            // Remove a and b from further consideration in this pass through locking.
            isLocked[max_A_index] = 1;
            isLocked[max_B_index] = 1;

            Lg[i_L] = gain;
            La[i_L] = a_index_A;
            Lb[i_L] = b_index_B;
            i_L++;

            // Update D values of affected elements in A - a and B - b
            cur = (graph->adjLists)[max_A_index];
            while (cur != NULL)
            {
                struct location loc = locations[cur->vertex];
                if (loc.a_or_b == 0)
                {
                    gain_A[loc.index] = gain_A[loc.index] + 1;
                }
                else
                {
                    gain_B[loc.index] = gain_B[loc.index] - 1;
                }
                next = cur->next;
                cur = next;
            }

            cur = (graph->adjLists)[max_B_index];
            while (cur != NULL)
            {
                struct location loc = locations[cur->vertex];
                if (loc.a_or_b == 1)
                {
                    gain_B[loc.index] = gain_B[loc.index] + 1;
                }
                else
                {
                    gain_A[loc.index] = gain_A[loc.index] - 1;
                }
                next = cur->next;
                cur = next;
            }
        }

        if (i_L == 0)
        {
            break;
        }

        // Find k which maximizes gmax
        int *gmaxes = malloc(size_of_A * sizeof(int));
        gmaxes[0] = Lg[0];
        for (i = 1; i < i_L; i++)
        {
            gmaxes[i] = gmaxes[i - 1] + Lg[i];
        }
        int k = 0;
        gmax = gmaxes[0];
        for (i = 1; i < i_L; i++)
        {
            if (gmaxes[i] > gmax)
            {
                gmax = gmaxes[i];
                k = i;
            }
        }

        int temp = 0;
        if (gmax > 0)
        {
            for (i = 0; i <= k; i++)
            {
                temp = A[La[i]];
                A[La[i]] = B[Lb[i]];
                B[Lb[i]] = temp;

                locations[B[Lb[i]]].a_or_b = 1;
                locations[A[La[i]]].a_or_b = 0;

                temp = locations[B[Lb[i]]].index;
                locations[B[Lb[i]]].index = locations[A[La[i]]].index;
                locations[A[La[i]]].index = temp;
            }
        }
        free(gmaxes);
    } while (gmax > 0);

    clock_t end = clock();
    double time_spent = 0.0;
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    *time_spentPtr = time_spent;

    // Calculate final cut value
    cut_value = calculate_cut_value(graph, A, size_of_A, B, size_of_B, locations);
    *final_cut_sizePtr = cut_value;

    // Clean
    free(gain_A);
    free(gain_B);
    free(Lg);
    free(La);
    free(Lb);

    free(isLocked);
    free(locations);
    free(A);
    free(B);

    return 1;
}

int simpler_KL_heap(struct Graph *graph, int *initial_cut_sizePtr, int *final_cut_sizePtr, double *time_spentPtr)
{
    // Split V into two balanced disjoint sets A and B
    int no_of_vertices = graph->numVertices;
    struct location *locations = malloc(no_of_vertices * sizeof(struct location));

    int size_of_A = no_of_vertices / 2;
    int *A = malloc(size_of_A * sizeof(int));
    int i;
    int j;

    for (i = 0; i < size_of_A; i++)
    {
        A[i] = i;
        locations[i].a_or_b = 0;
        locations[i].index = i;
    }

    int size_of_B = no_of_vertices - no_of_vertices / 2;
    int *B = malloc(size_of_B * sizeof(int));

    for (j = 0; j < size_of_B; j++)
    {
        B[j] = i;
        locations[i].a_or_b = 1;
        locations[i].index = j;
        i++;
    }

    // Calculate initial cut value
    int cut_value = calculate_cut_value(graph, A, size_of_A, B, size_of_B, locations);
    *initial_cut_sizePtr = cut_value;

    int gmax = 0;

    // Create lock array
    int *isLocked = malloc(no_of_vertices * sizeof(int));

    clock_t begin = clock();

    struct node *cur;
    struct node *next;

    // Loop mallocs
    struct heap_node *gain_A = malloc((size_of_A) * sizeof(struct heap_node));
    struct heap_node *gain_B = malloc((size_of_B) * sizeof(struct heap_node));

    int *heap_indexes_A = malloc(size_of_A * sizeof(int));
    int *heap_indexes_B = malloc(size_of_B * sizeof(int));

    int *Lg = malloc(size_of_A * sizeof(int));
    int *La = malloc(size_of_A * sizeof(int));
    int *Lb = malloc(size_of_A * sizeof(int));

    do
    {
        // Initilize lock array
        for (i = 0; i < no_of_vertices; i++)
            isLocked[i] = 0;

        // Compute D values (move gain values) for all A and B
        for (i = 0; i < size_of_A; i++)
        {
            gain_A[i].gain = 0;
            gain_A[i].index = i;
            heap_indexes_A[i] = i;
        }
        for (i = 0; i < size_of_B; i++)
        {
            gain_B[i].gain = 0;
            gain_B[i].index = i;
            heap_indexes_B[i] = i;
        }

        for (i = 0; i < size_of_A; i++)
        {
            cur = (graph->adjLists)[A[i]];
            while (cur != NULL)
            {
                struct location loc = locations[cur->vertex];

                if (loc.a_or_b == 0)
                {
                    gain_A[i].gain = gain_A[i].gain - 1;
                }
                else
                {
                    gain_A[i].gain = gain_A[i].gain + 1;
                }
                next = cur->next;
                cur = next;
            }
        }

        for (i = 0; i < size_of_B; i++)
        {
            cur = (graph->adjLists)[B[i]];
            while (cur != NULL)
            {
                struct location loc = locations[cur->vertex];

                if (loc.a_or_b == 1)
                {
                    gain_B[i].gain = gain_B[i].gain - 1;
                }
                else
                {
                    gain_B[i].gain = gain_B[i].gain + 1;
                }
                next = cur->next;
                cur = next;
            }
        }

        // Let Lg, La, Lb be empty lists and i = 0
        int i_L = 0;

        // Heapify the gain arrays
        int size_of_heap_A = size_of_A;
        int size_of_heap_B = size_of_B;

        build_heap(gain_A, heap_indexes_A, size_of_heap_A);
        build_heap(gain_B, heap_indexes_B, size_of_heap_B);

        int n;
        for (n = 0; n < size_of_A; n++)
        {
            // Find a and b such that the exchange gain g = Da + Db is maximized.
            if (size_of_heap_A <= 0 || size_of_B <= 0)
                break;

            struct heap_node max_gain_a = extract_max(gain_A, heap_indexes_A, &size_of_heap_A);
            struct heap_node max_gain_b = extract_max(gain_B, heap_indexes_B, &size_of_heap_B);

            while (isLocked[A[max_gain_a.index]] == 1)
            {
                if (size_of_heap_A == 0)
                    break;
                max_gain_a = extract_max(gain_A, heap_indexes_A, &size_of_heap_A);
            }

            if (isLocked[A[max_gain_a.index]] == 1)
                break;

            while (isLocked[B[max_gain_b.index]] == 1)
            {
                if (size_of_heap_B == 0)
                    break;
                max_gain_b = extract_max(gain_B, heap_indexes_B, &size_of_heap_B);
            }

            if (isLocked[B[max_gain_b.index]] == 1)
                break;

            int max_A_index = max_gain_a.index;
            int max_B_index = max_gain_b.index;

            int gain = max_gain_a.gain + max_gain_b.gain;
            int a_index_A = max_A_index;
            int b_index_B = max_B_index;

            max_A_index = A[max_A_index];
            max_B_index = B[max_B_index];

            cur = (graph->adjLists)[max_A_index];
            while (cur != NULL)
            {
                if (cur->vertex == max_B_index)
                {
                    gain = gain - 2;
                    break;
                }
                next = cur->next;
                cur = next;
            }

            // Remove a and b from further consideration in this pass through locking.
            isLocked[max_A_index] = 1;
            isLocked[max_B_index] = 1;

            Lg[i_L] = gain;
            La[i_L] = a_index_A;
            Lb[i_L] = b_index_B;
            i_L++;

            // Update D values of affected elements in A - a and B - b
            cur = (graph->adjLists)[max_A_index];
            while (cur != NULL)
            {
                struct location loc = locations[cur->vertex];
                if (loc.a_or_b == 0)
                {
                    int index_in_heap = heap_indexes_A[loc.index];
                    if (index_in_heap < size_of_heap_A)
                    {
                        int new_key = gain_A[index_in_heap].gain + 1;
                        increase_key(gain_A, heap_indexes_A, size_of_heap_A, index_in_heap, new_key);
                    }
                }
                else
                {
                    int index_in_heap = heap_indexes_B[loc.index];
                    if (index_in_heap < size_of_heap_B)
                    {
                        int new_key = gain_B[index_in_heap].gain - 1;
                        decrease_key(gain_B, heap_indexes_B, size_of_heap_B, index_in_heap, new_key);
                    }
                }
                next = cur->next;
                cur = next;
            }

            cur = (graph->adjLists)[max_B_index];
            while (cur != NULL)
            {
                struct location loc = locations[cur->vertex];
                if (loc.a_or_b == 1)
                {
                    int index_in_heap = heap_indexes_B[loc.index];
                    if (index_in_heap < size_of_heap_B)
                    {
                        int new_key = gain_B[index_in_heap].gain + 1;
                        increase_key(gain_B, heap_indexes_B, size_of_heap_B, index_in_heap, new_key);
                    }
                }
                else
                {
                    int index_in_heap = heap_indexes_A[loc.index];
                    if (index_in_heap < size_of_heap_A)
                    {
                        int new_key = gain_A[index_in_heap].gain - 1;
                        decrease_key(gain_A, heap_indexes_A, size_of_heap_A, index_in_heap, new_key);
                    }
                }
                next = cur->next;
                cur = next;
            }
        }

        if (i_L == 0)
        {
            break;
        }

        // Find k which maximizes gmax
        int *gmaxes = malloc(i_L * sizeof(int));
        gmaxes[0] = Lg[0];
        for (i = 1; i < i_L; i++)
        {
            gmaxes[i] = gmaxes[i - 1] + Lg[i];
        }
        int k = 0;
        gmax = gmaxes[0];
        for (i = 1; i < i_L; i++)
        {
            if (gmaxes[i] > gmax)
            {
                gmax = gmaxes[i];
                k = i;
            }
        }

        int temp = 0;
        if (gmax > 0)
        {
            for (i = 0; i <= k; i++)
            {
                temp = A[La[i]];
                A[La[i]] = B[Lb[i]];
                B[Lb[i]] = temp;

                locations[B[Lb[i]]].a_or_b = 1;
                locations[A[La[i]]].a_or_b = 0;

                temp = locations[B[Lb[i]]].index;
                locations[B[Lb[i]]].index = locations[A[La[i]]].index;
                locations[A[La[i]]].index = temp;
            }
        }
        free(gmaxes);
    } while (gmax > 0);

    clock_t end = clock();
    double time_spent = 0.0;
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    *time_spentPtr = time_spent;

    // Calculate final cut value
    cut_value = calculate_cut_value(graph, A, size_of_A, B, size_of_B, locations);
    *final_cut_sizePtr = cut_value;

    // Clean
    free(gain_A);
    free(gain_B);
    free(heap_indexes_A);
    free(heap_indexes_B);
    free(Lg);
    free(La);
    free(Lb);

    free(isLocked);
    free(locations);
    free(A);
    free(B);

    return 1;
}

// Main function
int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename] [a or b]\n", argv[0]);
        exit(1);
    }

    if (strcmp(argv[2], "a") != 0)
    {
        if (strcmp(argv[2], "b") != 0)
        {
            printf("%s\n", argv[2]);
            fprintf(stderr, "Usage: %s [martix-market-filename] [a or b]\n", argv[0]);
            exit(1);
        }
    }

    // Open and read file
    FILE *file;
    if ((file = fopen(argv[1], "r")) == NULL)
    {
        fprintf(stderr, "File %s cannot be opened\n", argv[1]);
        exit(1);
    }

    struct Graph *graph;
    int no_of_vertices = 0;

    if (readMatrixFile(file, &graph, &no_of_vertices) != 0)
    {
        fprintf(stderr, "Could not read file %s\n", argv[1]);
        exit(1);
    }

    if (strcmp(argv[2], "a") == 0)
    {
        // Simpler KL runs
        // Heap implementation
        int initial_cut_size_heap = 0;
        int final_cut_size_heap = 0;
        double time_spent_heap = 0.0;

        simpler_KL_heap(graph, &initial_cut_size_heap, &final_cut_size_heap, &time_spent_heap);
        printf("%d %d %fs\n", initial_cut_size_heap, final_cut_size_heap, time_spent_heap);
    }
    else
    {
        // Straightforward implementation
        int initial_cut_size = 0;
        int final_cut_size = 0;
        double time_spent = 0.0;

        simpler_KL(graph, &initial_cut_size, &final_cut_size, &time_spent);
        printf("%d %d %fs\n", initial_cut_size, final_cut_size, time_spent);
    }

    int i;
    struct node **lists = graph->adjLists;
    for (i = 0; i < no_of_vertices; i++)
    {
        struct node *cur = lists[i];
        struct node *next;
        while (cur != NULL)
        {
            next = cur->next;
            cur->next = NULL;
            free(cur);
            cur = next;
        }
    }
    return 0;
}
