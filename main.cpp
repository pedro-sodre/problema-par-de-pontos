// Franco Barpp Gomes - 2126613
// Pedro Sodré dos Santos - 2126745

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <iomanip>
#include <limits>
#include <chrono>
#include <cmath>

// #define TESTING

const double INFTY = std::numeric_limits<double>::infinity();

struct Coord {
    double x;
    double y;
};

using CoordPair = std::pair<Coord, Coord>;

std::vector<Coord> parseFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    std::string lineBuffer;
    std::vector<Coord> coordinates;

    std::getline(inputFile, lineBuffer);
    const std::size_t length = std::stoi(lineBuffer);

    coordinates.reserve(length);

    while (std::getline(inputFile, lineBuffer)) {
        std::istringstream lineStream(lineBuffer);
        Coord coord;

        lineStream >> std::setprecision(6) >> std::fixed >> coord.x >> coord.y;
        coordinates.push_back(coord);
    }
    return coordinates;
}

double euclideanDistance(CoordPair cp) {
    const Coord c1 = cp.first;
    const Coord c2 = cp.second;
    return sqrt((c1.x - c2.x)*(c1.x - c2.x) + (c1.y - c2.y)*(c1.y - c2.y));
}

CoordPair returnClosestStrip(std::vector<Coord> coordinates, std::size_t from, std::size_t to, Coord splitPos, double delta) {
    std::vector<Coord> nearStrip;

    std::copy_if(coordinates.begin() + from, coordinates.begin() + to + 1, std::back_inserter(nearStrip), [splitPos, delta](Coord c) {
        return std::abs(c.x - splitPos.x) < delta;
    });

    std::sort(nearStrip.begin(), nearStrip.end(), [](Coord c1, Coord c2) {
        return c1.y < c2.y;
    });

    CoordPair bestPair({INFTY, INFTY}, {-INFTY, -INFTY});
    double bestDistance = INFTY;

    for (std::size_t i = 0; i < nearStrip.size(); i++) {
        for (std::size_t j = i + 1; j < nearStrip.size(); j++) {
            Coord c1 = nearStrip[i];
            Coord c2 = nearStrip[j];
            CoordPair cp = CoordPair(c1, c2);
            
            if ((c2.y - c1.y) > bestDistance) {
                break;
            }
            
            const double dist = euclideanDistance(cp);
            
            if (dist < bestDistance) {
                bestDistance = dist;
                bestPair = cp;
            }
        }
    }
    
    return bestPair;
}

CoordPair returnClosest(std::vector<Coord> coordinates, std::size_t from, std::size_t to) {
    const std::size_t size = to - from + 1;

    if (size > 3) {
        const CoordPair left = returnClosest(coordinates, from, from + size / 2);
        const CoordPair right = returnClosest(coordinates, from + size / 2, to);

        const double d1 = euclideanDistance(left);
        const double d2 = euclideanDistance(right);

        const CoordPair center = returnClosestStrip(coordinates, from, to, coordinates[from + size / 2], std::min(d1, d2));
        const double d3 = euclideanDistance(center);

        if (d1 < d2) return (d1 < d3)? left : center;
        else return (d2 < d3)? right : center;
    }
    else {
        switch (size) {
            case 1:
                throw std::runtime_error("No pair!");
                break;
            case 2:
                return CoordPair(coordinates[from], coordinates[to]);
                break;
            case 3: {
                return std::min({
                    CoordPair(coordinates[from], coordinates[from + 1]),
                    CoordPair(coordinates[from], coordinates[from + 2]),
                    CoordPair(coordinates[from + 1], coordinates[from + 2]),
                }, [](CoordPair p1, CoordPair p2) {
                    return euclideanDistance(p1) < euclideanDistance(p2);
                });
                break;
            }
        }
    }

    return {{0, 0}, {0, 0}};
}

CoordPair closestPairs(std::vector<Coord> coords) {
    std::sort(coords.begin(), coords.end(), [](Coord c1, Coord c2) {
        return c1.x < c2.x;
    });
    return returnClosest(coords, 0, coords.size() - 1);
}

// Algoritmo O(n^2) para teste
CoordPair testClosestPairs(std::vector<Coord> coords) {
    double bestDistance = std::numeric_limits<double>::max();
    CoordPair bestPair;
    
    for (std::size_t i = 0; i < coords.size(); i++) {
        for (std::size_t j = i + 1; j < coords.size(); j++) {
            const CoordPair cp(coords[i], coords[j]);
            const double dist = euclideanDistance(cp);

            if (dist < bestDistance) {
                bestDistance = dist;
                bestPair = cp;
            }
        }
    }

    return bestPair;
}

int main(int argc, char** argv) {
    if (argc != 2) return 1;

    std::string filename = argv[1];
    std::vector<Coord> list = parseFile(filename);

    const auto now = std::chrono::high_resolution_clock::now();
    CoordPair best = closestPairs(list);
    double dist = euclideanDistance(best);
    
    const auto after = std::chrono::high_resolution_clock::now();
    const double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(after - now).count() / 1e9;

    std::cout << std::setprecision(6) << std::fixed << duration <<  " " << dist << " " << best.first.x << " " << best.first.y << " " << best.second.x << " " << best.second.y << std::endl;

#ifdef TESTING
    std::vector<Coord> list2 = parseFile(filename);
    CoordPair best2 = testClosestPairs(list2);
    double dist2 = euclideanDistance(best2);

    std::cout << "Test: " << std::setprecision(6) << std::fixed << duration <<  " " << dist2 << " " << best2.first.x << " " << best2.first.y << " " << best2.second.x << " " << best2.second.y << std::endl;

    if (dist != dist2) {
        throw std::runtime_error("Results should be the same");
    }
#endif
    
    return 0;
}

