//
// Created by ssunah on 6/22/18.
//

#include "graph.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <utility/graphoperations.h>
#include <utility/graphoperations.cpp>

// Function to build a reverse index for the graph.
void Graph::BuildReverseIndex() {
    // Allocate memory for the reverse index array based on the number of vertices.
    reverse_index_ = new ui[vertices_count_];
    // Allocate memory for the offsets array based on the number of labels plus one.
    reverse_index_offsets_= new ui[labels_count_ + 1];
    // Initialize the first offset to 0.
    reverse_index_offsets_[0] = 0;

    // Variable to keep track of the total number of vertices processed.
    ui total = 0;
    // Loop over each label to calculate offsets.
    for (ui i = 0; i < labels_count_; ++i) {
        // Set the offset for the current label to the total number of vertices processed so far.
        reverse_index_offsets_[i + 1] = total;
        // Increment the total by the frequency of the current label.
        total += labels_frequency_[i];
    }

    // Loop over all vertices to build the reverse index.
    for (ui i = 0; i < vertices_count_; ++i) {
        // Get the label of the current vertex.
        LabelID label = labels_[i];
        // Place the vertex at the next available position for its label in the reverse index.
        // Then, increment the offset for the label to point to the next available position.
        reverse_index_[reverse_index_offsets_[label + 1]++] = i;
    }
}




// Conditionally compile this function only if OPTIMIZED_LABELED_GRAPH is set to 1.
#if OPTIMIZED_LABELED_GRAPH == 1

// Function to build the Neighbor Label Frequency (NLF) structure for the graph.
void Graph::BuildNLF() {
    // Allocate memory for the NLF array, which is an array of unordered maps.
    // Each map will hold the frequency of each label among the neighbors of a vertex.
    nlf_ = new std::unordered_map<LabelID, ui>[vertices_count_];

    // Iterate over all vertices in the graph.
    for (ui i = 0; i < vertices_count_; ++i) {
        // Variable to store the number of neighbors.
        ui count;
        // Get the neighbors of the current vertex.
        const VertexID * neighbors = getVertexNeighbors(i, count);

        // Iterate over each neighbor of the current vertex.
        for (ui j = 0; j < count; ++j) {
            // Get the ID of the current neighbor.
            VertexID u = neighbors[j];
            // Get the label of the current neighbor.
            LabelID label = getVertexLabel(u);
            // If the label is not already in the map for the current vertex, initialize its frequency to 0.
            if (nlf_[i].find(label) == nlf_[i].end()) {
                nlf_[i][label] = 0;
            }

            // Increment the frequency of the label for the current vertex.
            nlf_[i][label] += 1;
        }
    }
}


// Function to build the label offset structure for the graph.
void Graph::BuildLabelOffset() {
    // Calculate the size needed for the labels_offsets_ array.
    size_t labels_offset_size = (size_t)vertices_count_ * labels_count_ + 1;
    // Allocate memory for the labels_offsets_ array.
    labels_offsets_ = new ui[labels_offset_size];
    // Initialize the labels_offsets_ array with zeros.
    std::fill(labels_offsets_, labels_offsets_ + labels_offset_size, 0);

    // Sort the neighbors of each vertex by their labels (and by their IDs if labels are equal).
    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1],
            [this](const VertexID u, const VertexID v) -> bool {
                return labels_[u] == labels_[v] ? u < v : labels_[u] < labels_[v];
            });
    }

    // Build the labels_offsets_ array with the correct offsets.
    for (ui i = 0; i < vertices_count_; ++i) {
        // Initialize the previous and current label variables.
        LabelID previous_label = 0;
        LabelID current_label = 0;

        // Calculate the starting index in the labels_offsets_ array for the current vertex.
        labels_offset_size = i * labels_count_;
        // Set the initial offset for the current vertex.
        labels_offsets_[labels_offset_size] = offsets_[i];

        // Iterate over the neighbors of the current vertex.
        for (ui j = offsets_[i]; j < offsets_[i + 1]; ++j) {
            // Get the label of the current neighbor.
            current_label = labels_[neighbors_[j]];

            // If the label changes, update the offsets for the range of labels.
            if (current_label != previous_label) {
                for (ui k = previous_label + 1; k <= current_label; ++k) {
                    labels_offsets_[labels_offset_size + k] = j;
                }
                // Update the previous label to the current label.
                previous_label = current_label;
            }
        }

        // After processing all neighbors, set the offsets for the remaining labels.
        for (ui l = current_label + 1; l <= labels_count_; ++l) {
            labels_offsets_[labels_offset_size + l] = offsets_[i + 1];
        }
    }
}

#endif // This ends the conditional compilation block.


// Function to load a graph from a file.
void Graph::loadGraphFromFile(const std::string &file_path) {
    // Open the file at the given file path.
    std::ifstream infile(file_path);

    // Check if the file is successfully opened.
    if (!infile.is_open()) {
        std::cout << "Cannot open the graph file " << file_path << " ." << std::endl;
        // Exit if the file cannot be opened.
        exit(-1);
    }

    // Variables to store the type of the line being read, and the counts of vertices and edges.
    char type;
    infile >> type >> vertices_count_ >> edges_count_;
    // Allocate memory for the offsets array and initialize the first element to 0.
    offsets_ = new ui[vertices_count_ +  1];
    offsets_[0] = 0;

    // Allocate memory for the neighbors array and the labels array.
    neighbors_ = new VertexID[edges_count_ * 2];
    labels_ = new LabelID[vertices_count_];
    // Initialize the labels count and the maximum degree of the graph.
    labels_count_ = 0;
    max_degree_ = 0;

    // Variable to keep track of the maximum label ID.
    LabelID max_label_id = 0;
    // Vector to keep track of the current offset for each vertex's neighbors.
    std::vector<ui> neighbors_offset(vertices_count_, 0);

    // Read the file line by line.
    while (infile >> type) {
        if (type == 'v') { // If the line defines a vertex.
            VertexID id;
            LabelID label;
            ui degree;
            infile >> id >> label >> degree;

            // Set the label for the vertex and update the offsets array.
            labels_[id] = label;
            offsets_[id + 1] = offsets_[id] + degree;

            // Update the maximum degree if necessary.
            if (degree > max_degree_) {
                max_degree_ = degree;
            }

            // Update the labels frequency map.
            if (labels_frequency_.find(label) == labels_frequency_.end()) {
                labels_frequency_[label] = 0;
                // Update the maximum label ID if necessary.
                if (label > max_label_id)
                    max_label_id = label;
            }

            // Increment the frequency of the label.
            labels_frequency_[label] += 1;
        }
        else if (type == 'e') { // If the line defines an edge.
            VertexID begin;
            VertexID end;
            infile >> begin >> end;

            // Add the edge to the neighbors array for both vertices.
            ui offset = offsets_[begin] + neighbors_offset[begin];
            neighbors_[offset] = end;

            offset = offsets_[end] + neighbors_offset[end];
            neighbors_[offset] = begin;

            // Increment the neighbors offset for both vertices.
            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;
        }
    }

    // Close the file after reading.
    infile.close();
    // Determine the labels count based on the maximum label ID and the size of the labels frequency map.
    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    // Find the maximum label frequency.
    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    // Sort the neighbors of each vertex.
    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }

    // Build the reverse index for the graph.
    BuildReverseIndex();

    // If the graph is optimized for labeled operations and the label offset is enabled, build the NLF.
#if OPTIMIZED_LABELED_GRAPH == 1
    if (enable_label_offset_) {
        BuildNLF();
        // The BuildLabelOffset function is commented out, indicating it may be optional or not needed.
        // BuildLabelOffset();
    }
#endif
}


void Graph::printGraphMetaData() {
    std::cout << "|V|: " << vertices_count_ << ", |E|: " << edges_count_ << ", |\u03A3|: " << labels_count_ << std::endl;
    std::cout << "Max Degree: " << max_degree_ << ", Max Label Frequency: " << max_label_frequency_ << std::endl;
}

void Graph::buildCoreTable() {
    core_table_ = new int[vertices_count_];
    GraphOperations::getKCore(this, core_table_);

    for (ui i = 0; i < vertices_count_; ++i) {
        if (core_table_[i] > 1) {
            core_length_ += 1;
        }
    }
}

void Graph::loadGraphFromFileCompressed(const std::string &degree_path, const std::string &edge_path,
                                        const std::string &label_path) {
    std::ifstream deg_file(degree_path, std::ios::binary);

    if (deg_file.is_open()) {
        std::cout << "Open degree file " << degree_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot open degree file " << degree_path << " ." << std::endl;
        exit(-1);
    }

    auto start = std::chrono::high_resolution_clock::now();
    int int_size;
    deg_file.read(reinterpret_cast<char *>(&int_size), 4);
    deg_file.read(reinterpret_cast<char *>(&vertices_count_), 4);
    deg_file.read(reinterpret_cast<char *>(&edges_count_), 4);

    offsets_ = new ui[vertices_count_ + 1];
    ui* degrees = new unsigned int[vertices_count_];

    deg_file.read(reinterpret_cast<char *>(degrees), sizeof(int) * vertices_count_);


    deg_file.close();
    deg_file.clear();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Load degree file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;

    std::ifstream adj_file(edge_path, std::ios::binary);

    if (adj_file.is_open()) {
        std::cout << "Open edge file " << edge_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot open edge file " << edge_path << " ." << std::endl;
        exit(-1);
    }

    start = std::chrono::high_resolution_clock::now();
    size_t neighbors_count = (size_t)edges_count_ * 2;
    neighbors_ = new ui[neighbors_count];

    offsets_[0] = 0;
    for (ui i = 1; i <= vertices_count_; ++i) {
        offsets_[i] = offsets_[i - 1] + degrees[i - 1];
    }

    max_degree_ = 0;

    for (ui i = 0; i < vertices_count_; ++i) {
        if (degrees[i] > 0) {
            if (degrees[i] > max_degree_)
                max_degree_ = degrees[i];
            adj_file.read(reinterpret_cast<char *>(neighbors_ + offsets_[i]), degrees[i] * sizeof(int));
            std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
        }
    }

    adj_file.close();
    adj_file.clear();

    delete[] degrees;

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Load adj file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;


    std::ifstream label_file(label_path, std::ios::binary);
    if (label_file.is_open())  {
        std::cout << "Open label file " << label_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot open label file " << label_path << " ." << std::endl;
        exit(-1);
    }

    start = std::chrono::high_resolution_clock::now();

    labels_ = new ui[vertices_count_];
    label_file.read(reinterpret_cast<char *>(labels_), sizeof(int) * vertices_count_);

    label_file.close();
    label_file.clear();

    ui max_label_id = 0;
    for (ui i = 0; i < vertices_count_; ++i) {
        ui label = labels_[i];

        if (labels_frequency_.find(label) == labels_frequency_.end()) {
            labels_frequency_[label] = 0;
            if (label > max_label_id)
                max_label_id = label;
        }

        labels_frequency_[label] += 1;
    }

    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Load label file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    BuildReverseIndex();
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Build reverse index file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
#if OPTIMIZED_LABELED_GRAPH == 1
    if (enable_label_offset_) {
        BuildNLF();
        // BuildLabelOffset();
    }
#endif
}

void Graph::storeComparessedGraph(const std::string& degree_path, const std::string& edge_path,
                                  const std::string& label_path) {
    ui* degrees = new ui[vertices_count_];
    for (ui i = 0; i < vertices_count_; ++i) {
        degrees[i] = offsets_[i + 1] - offsets_[i];
    }

    std::ofstream deg_outputfile(degree_path, std::ios::binary);

    if (deg_outputfile.is_open()) {
        std::cout << "Open degree file " << degree_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot degree edge file " << degree_path << " ." << std::endl;
        exit(-1);
    }

    int int_size = sizeof(int);
    size_t vertex_array_bytes = ((size_t)vertices_count_) * 4;
    deg_outputfile.write(reinterpret_cast<const char *>(&int_size), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(&vertices_count_), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(&edges_count_), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(degrees), vertex_array_bytes);

    deg_outputfile.close();
    deg_outputfile.clear();

    delete[] degrees;

    std::ofstream edge_outputfile(edge_path, std::ios::binary);

    if (edge_outputfile.is_open()) {
        std::cout << "Open edge file " << edge_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot edge file " << edge_path << " ." << std::endl;
        exit(-1);
    }

    size_t edge_array_bytes = ((size_t)edges_count_ * 2) * 4;
    edge_outputfile.write(reinterpret_cast<const char *>(neighbors_), edge_array_bytes);

    edge_outputfile.close();
    edge_outputfile.clear();

    std::ofstream label_outputfile(label_path, std::ios::binary);

    if (label_outputfile.is_open()) {
        std::cout << "Open label file " << label_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot label file " << label_path << " ." << std::endl;
        exit(-1);
    }

    size_t label_array_bytes = ((size_t)vertices_count_) * 4;
    label_outputfile.write(reinterpret_cast<const char *>(labels_), label_array_bytes);

    label_outputfile.close();
    label_outputfile.clear();
}
