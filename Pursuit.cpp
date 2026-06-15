#include <cmath>
#include <stack>
#include <numbers>
#include <algorithm>
#include <iostream>

#include <nlohmann/json.hpp>

#include "Leonard.hpp"
#include "Boundary.hpp"
#include "Modal.hpp"

using json = nlohmann::json;

Qualify_Results Leonard::do_Qualify(const std::vector<float> *data, const Measurement X, const std::vector<bool> &included, std::string pop_name)
{
    Qualify_Results results(X);
    thread_local std::vector<float> x;
    x.clear();
    x.reserve(included.size());

    // get statistics for this measurement for this subset
    double Sx = 0, Sxx = 0;
    size_t n = 0, m = 0;
    auto it = included.begin();
    while ((it = std::find(it, included.end(), true)) != included.end())
    {
        size_t i = it - included.begin();
        float value = (*data)[i];
        ++n;
        Sx += value;
        Sxx += value * value;
        x.push_back(value);
        ++it;
    }
    std::sort(x.begin(), x.end());
    x[n] = 1;
    while (x[m] == 0)
        ++m; // for CyToF/exponential we censor true zeros
    if (m == n)
        return results;
    const double mu = Sx / (double)n;
    const double sigma = sqrt((Sxx - Sx * mu) / (double)(n - 1));
    const double lambda = Sx / (double)(n - m);

    // compute Kullback-Leibler Divergence
    if (sigma > 0)
    {
        // normalization factors for truncated distribution
        double NQn = .5 * (erf((x[n] - mu) / sigma / std::numbers::sqrt2) - erf((x[0] - mu) / sigma / std::numbers::sqrt2));
        double NQe = exp(-x[0] / lambda) - exp(-x[n] / lambda);
        for (size_t i = 0, j; i < n; i = j)
        {
            j = i + 1;
            while ((x[j] - x[i]) < .001 && j < n)
                j++;

            double P = (double)(j - i) / (double)n;
            double Q = .5 * (erf((x[j] - mu) / sigma / std::numbers::sqrt2) - erf((x[i] - mu) / sigma / std::numbers::sqrt2)) / NQn;
            if (Q > 0)                  // catch underflow that causes infinite result
                results.KLDn += P * log(P / Q); // I didn't think it was possible either

            if (i == 0 && m > 0)
                P = (double)(j - m) / (double)(n - m);
            else
                P = (double)(j - i) / (double)(n - m);
            Q = (exp(-x[i] / mu) - exp(-x[j] / mu)) / NQe;
            if (Q > 0)
                results.KLDe += P * log(P / Q);
        }
    }
    results.qualified = results.KLDn > selections.kld_norm && results.KLDe > selections.kld_exp;
    return results;
}

Projection_Results Leonard::do_Projection(const std::vector<std::vector<float> *> &data, const Measurement X, const Measurement Y, const std::vector<bool> &included, std::string pop_name)
{
    Projection_Results candidate(X,Y);
    Weighty<2> weighty(256);
    Event<2> event;
    for (auto it = std::find(included.begin(), included.end(), true); it != included.end(); it = std::find(++it, included.end(), true))
    {
        event[0] = (*data[candidate.X])[it - included.begin()];
        event[1] = (*data[candidate.Y])[it - included.begin()];
        weighty.event(event);
    }

    std::unique_ptr<Kernel<2, float>> var_kernel = std::make_unique<Kernel<2, float>>(weighty);
    ModalClustering modal;
    ColoredBoundary cluster_bounds;
    std::vector<ColoredEdge> edges;
    do
    {
        do
        {
            if (candidate.pass == 10)
            {
                candidate.outcome = Projection_Results::Status::EPP_error;
                return candidate;
            }
            // last density becomes this variance estimator
            std::swap(weighty.kernel, var_kernel);
            weighty.kernel->radius(selections.smoothing * std::pow(std::numbers::sqrt2, ++candidate.pass));
            // apply kernel to cosine transform
            weighty.transform(weighty.weight, weighty.cosine);
            weighty.apply_kernel(weighty.cosine, weighty.filtered, weighty.kernel.get());
            // inverse discrete cosine transform
            // gives a smoothed density estimator
            weighty.transform(weighty.filtered, weighty.density);
            // modal clustering
            candidate.clusters = modal.findClusters(weighty.density.data, candidate.pass, selections);
        } while (candidate.clusters > 12);
        if (candidate.clusters < 2)
        {
            candidate.outcome = Projection_Results::Status::EPP_no_cluster;
            return candidate;
        }

        modal.getBoundary(weighty.density.data, cluster_bounds);
        // get the edges, which have their own weights
        edges = cluster_bounds.getEdges();
        // smooth some more if graph is too complex to process
    } while (edges.size() > max_booleans);

        // get the dual graph of the map
    ColoredGraph graph = cluster_bounds.getDualGraph();

    // Density Based Merging
    Function<2, float> variance(weighty);
    if (candidate.pass == 1)
    { // otherwise it was swapped in above
        var_kernel->radius(selections.smoothing);
        weighty.apply_kernel(weighty.cosine, weighty.filtered, var_kernel.get());
        weighty.transform(weighty.filtered, variance);
    }
    thread_local struct cluster_merge
    {
        float edge_max, edge_var;
        BitPosition i;

        bool operator<(
            const cluster_merge &other) const noexcept
        { // larger taken first so sense inverted
            return edge_max > other.edge_max;
        };
    } merge[max_booleans];
    // find the highest point on each edge
    for (BitPosition i = 0; i < edges.size(); ++i)
    {
        const ColoredEdge &edge = edges[i];
        ColoredPoint point = edge.points[0];
        float edge_max = weighty.density[point.i + (N + 1) * point.j];
        for (Count j = 1; j < edge.points.size(); ++j)
        {
            ColoredPoint p = edge.points[j];
            float d = weighty.density[p.i + (N + 1) * p.j];
            if (d > edge_max)
            {
                point = p;
                edge_max = d;
            }
        }
        merge[i].edge_max = edge_max;
        merge[i].edge_var = variance[point.i + (N + 1) * point.j];
        merge[i].i = i;
    }
    std::sort(merge, merge + edges.size());
    // now test them in decreasong order of height
    for (BitPosition j = 0; j < edges.size(); ++j)
    {
        const ColoredEdge &edge = edges[merge[j].i];
        // can't apply it to half edges of border square
        if (edge.clockwise == 0 || edge.widdershins == 0)
            continue;

        double cluster_max, cluster_var;
        // the smaller of the maxima of the clusters the edge divides
        if (modal.maxima[edge.clockwise] < modal.maxima[edge.widdershins])
        {
            cluster_max = modal.maxima[edge.clockwise];
            ColoredPoint point = modal.center[edge.clockwise];
            cluster_var = variance[point.i + (N + 1) * point.j];
        }
        else
        {
            cluster_max = modal.maxima[edge.widdershins];
            ColoredPoint point = modal.center[edge.widdershins];
            cluster_var = variance[point.i + (N + 1) * point.j];
        }
        // formulas from DBM paper. 4N^2 normalizes the FFT
        // phi^2 is also gausian with sqrt(2) narrower kernel, but doesn't integrate to 1
        // compute with the normalized kernel and then 4 * radius^2 * pi corrects the integral
        double radius = N * selections.smoothing * std::pow(std::numbers::sqrt2, candidate.pass);
        double f_e = merge[j].edge_max / 4 / N / N / weighty.size();
        double sigma2_e = (merge[j].edge_var / 4 / N / N / 4 / radius / radius / pi / weighty.size() - f_e * f_e) / (weighty.size() - 1);
        double f_c = cluster_max / 4 / N / N / weighty.size();
        double sigma2_c = (cluster_var / 4 / N / N / 4 / radius / radius / std::numbers::pi / weighty.size() - f_c * f_c) / (weighty.size() - 1);
        // if the dip isn't significant, merge the two clusters and remove the edge
        if (sigma2_e > 0 && sigma2_c > 0 && f_c - std::sqrt(sigma2_c) < f_e + std::sqrt(sigma2_e))
        {
            // this may cause a lower one to become more significant
            if (modal.maxima[edge.clockwise] < modal.maxima[edge.widdershins])
            {
                modal.maxima[edge.clockwise] = modal.maxima[edge.widdershins];
                modal.center[edge.clockwise] = modal.center[edge.widdershins];
            }
            else
            {
                modal.maxima[edge.widdershins] = modal.maxima[edge.clockwise];
                modal.center[edge.widdershins] = modal.center[edge.clockwise];
            }

            graph = graph.simplify(merge[j].i);
            ++candidate.merges;
        }
    }
    // make sure there's anything left
    if (graph.isTrivial())
    {
        candidate.outcome = Projection_Results::Status::EPP_not_significant;
        return candidate;
    }

    // compute the cluster weights
    auto cluster_map = cluster_bounds.getMap();

    // find and score simple sub graphs
    struct
    {
        double score = std::numeric_limits<double>::infinity();
        Booleans edges;
        Booleans clusters;
        double edge_weight;
    } best;
    // pile of graphs to consider
    std::stack<ColoredGraph> pile;
    pile.push(graph);
    while (!pile.empty())
    {
        ++candidate.graphs;
        ColoredGraph graph = pile.top();
        pile.pop();
        if (graph.isSimple()) // one edge, i.e., two populations
        {
            // because the mode is always in the first cluster
            // we choose the node that includes it, i.e.,
            // the "in" set will always include the sample mode
            Booleans in_clusters = graph.left() & 1 ? graph.left() : graph.right();

            Booleans dual_edges = graph.edge();
            double edge_weight = 0;
            for (size_t i = 0; i < edges.size(); i++)
            {
                if (dual_edges & (1 << i))
                    edge_weight += edges[i].weight;
            }
            edge_weight /= 8 * N * N; // approximates number of events within a border region of width W;
            double score = edge_weight;
            double balance_factor = 0;
            unsigned long int in_weight = 0;
            assert(score > 0);
            // score this separatrix
            if (score < best.score)
            {
                best.score = score;
                best.edges = dual_edges;
                best.clusters = in_clusters;
                best.edge_weight = edge_weight;
            }
        }
        else
        { // not simple so simplify it some, i.e., remove one dual edge at a time
            // and merge two adjacent subsets. that makes a bunch more graphs to look at
            std::vector<ColoredGraph> simplified = graph.simplify();
            for (const auto &graph : simplified)
                pile.push(graph);
        }
    }
    if (best.score == std::numeric_limits<double>::infinity())
    {
        candidate.outcome = Projection_Results::Status::EPP_no_cluster;
        return candidate;
    }

    // find the separatrix
    ColoredBoundary subset_boundary;
    subset_boundary.clear();
    std::vector<ColoredPoint> interior_vertex;
    assert(interior_vertex.empty());
    std::vector<BitPosition> half_edges;
    half_edges.clear();
    for (BitPosition i = 0; i < edges.size(); i++)
    {
        ColoredEdge &edge = edges[i];
        // skip half edges from vertex squares for now
        if (edge.clockwise == 0 || edge.widdershins == 0)
        {
            half_edges.push_back(i);
            continue;
        }
        if (best.edges & (1 << i))
        {
            bool lefty = best.clusters & (1 << (edge.widdershins - 1));
            subset_boundary.addEdge(edge.points, !lefty, lefty);
            // end points on the boundaries of data space are vertices
            ColoredPoint point = edge.points.front();
            if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                subset_boundary.addVertex(point);
            else
            {
                auto position = std::find(interior_vertex.begin(), interior_vertex.end(), point);
                if (position != interior_vertex.end())
                    interior_vertex.erase(position);
                else
                    interior_vertex.push_back(point);
            }
            point = edge.points.back();
            if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                subset_boundary.addVertex(point);
            else
            {
                auto position = std::find(interior_vertex.begin(), interior_vertex.end(), point);
                if (position != interior_vertex.end())
                    interior_vertex.erase(position);
                else
                    interior_vertex.push_back(point);
            }
        }
    }
    while (!interior_vertex.empty())
    {   // need to fill in using some of the half edges
        bool making_progress = false;
        ColoredPoint end_point = interior_vertex.back(); 
        // so it will be the last one pushed bellow if any
        // there are two half edges and we must always extend the
        // same one to completion or it can deadlock or switchback
        for (auto hep = half_edges.begin(); hep != half_edges.end(); ++hep)
        {
            ColoredEdge &edge = edges[*hep];
            if (edge.points.front() != end_point && edge.points.back() != end_point)
                continue;
            // found the next connecting piece
            bool lefty;
            // since there are now only two, we can compute the missing color
            if (edge.widdershins == 0)
                lefty = !(best.clusters & (1 << (edge.clockwise - 1)));
            else 
                lefty = best.clusters & (1 << (edge.widdershins - 1));
            subset_boundary.addEdge(edge.points, !lefty, lefty);
            ColoredPoint point = edge.points.front();
            if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                subset_boundary.addVertex(point);
            else
            {
                auto position = std::find(interior_vertex.begin(), interior_vertex.end(), point);
                if (position != interior_vertex.end())
                    interior_vertex.erase(position);
                else
                    interior_vertex.push_back(point);
            }
            point = edge.points.back();
            if (point.i == 0 || point.i == N || point.j == 0 || point.j == N)
                subset_boundary.addVertex(point);
            else
            {
                auto position = std::find(interior_vertex.begin(), interior_vertex.end(), point);
                if (position != interior_vertex.end())
                    interior_vertex.erase(position);
                else
                    interior_vertex.push_back(point);
            }
            // don't use this one again
            half_edges.erase(hep);
            making_progress = true;
            break;
        }
        assert(making_progress);
    }
    subset_boundary.setColorful(2);
    assert(!subset_boundary.empty());

    ColoredEdge separatrix = subset_boundary.getEdges().front();
    candidate.separatrix.reserve(separatrix.points.size());
    for (ColoredPoint cp : separatrix.points)
        candidate.separatrix.push_back(cp);
    if (separatrix.widdershins)
        std::reverse(candidate.separatrix.begin(), candidate.separatrix.end());

    // create in/out subsets
    auto subset_map = subset_boundary.getMap();
    candidate.in.set.resize(weighty.size());
    candidate.out.set.resize(weighty.size());
    std::fill(candidate.in.set.begin(), candidate.in.set.end(), false);
    std::fill(candidate.out.set.begin(), candidate.out.set.end(), false);

    for (auto it = std::find(included.begin(), included.end(), true); it != included.end(); it = std::find(++it, included.end(), true))
        {
            size_t event = it - included.begin();
            
            double x = (*data[candidate.X])[event];
            double y = (*data[candidate.Y])[event];
            bool member = subset_map->colorAt(x, y);
            if (member)
            {
                ++candidate.in.count;
                candidate.in.set[event] = true;
            }
            else
            {
                ++candidate.out.count;
                candidate.out.set[event] = true;
            }
        }
    // and the check is in the mail
    candidate.outcome = Projection_Results::Status::EPP_success;
    candidate.score = best.score;

    return candidate;
}

Pursuit_Results Leonard::do_Pursuit(const std::vector<std::vector<float> *> &data, const std::vector<bool> &included, std::string pop_name)
{
    Pursuit_Results results;
    std::vector<std::unique_ptr<Qualify_Results>> qualifications;
    std::vector<std::future<Qualify_Results>> future_qualify;
    qualifications.reserve(data.size());
    future_qualify.reserve(data.size());
    for (size_t i = 0; i < data.size(); ++i)
        future_qualify.push_back(compute_plane.enqueue([this, i, data, included, pop_name]() {
                    return do_Qualify(data[i], i, included, pop_name);}));
    for (auto &result : future_qualify)
    {
        Qualify_Results res = result.get();
        qualifications.push_back(std::make_unique<Qualify_Results>(res));
    }
    std::sort(qualifications.begin(), qualifications.end(), [](const std::unique_ptr<Qualify_Results> &a, const std::unique_ptr<Qualify_Results> &b) 
        {   return a->KLDn > b->KLDn; });
    for (auto &q : qualifications)
        if (q->qualified)
            results.qualified.push_back(q->X);
    std::vector<std::future<Projection_Results>> future_projection;
    if (results.qualified.size() < 2)
    {
        results.qualified.push_back(qualifications[0]->X);
        results.qualified.push_back(qualifications[1]->X);
    }
    else
        for (unsigned i = 1; i < results.qualified.size(); ++i)
            for (unsigned j = 0; j < i; ++j)
            {
                std::vector<std::vector<float> *> plane = {data[j], data[i]};
                future_projection.push_back(compute_plane.enqueue([this, plane, i, j, included, pop_name]() {
                        return do_Projection(plane, j, i, included, pop_name); }));
            }
    for (auto &result : future_projection)
    {
        auto res = std::make_unique<Projection_Results>(result.get());
        if (!results.best_split || res->score < results.best_score)
        {
            results.best_score = res->score;
            results.best_split = std::move(res);
        }
    }
    if (results.best_split)
    {
        size_t min_count = static_cast<size_t>(selections.min_cluster_rel * included.size());
        if (min_count < selections.min_events)
            min_count = selections.min_events;
        if (results.best_split->outcome == Projection_Results::Status::EPP_success)
        {
            if (results.best_split->in.count >= min_count)
                results.future_children.push_back(control_plane.enqueue([this, data, in_set = results.best_split->in.set, pop_name]() {
                                return do_Pursuit(data, in_set, pop_name); }));
            if (results.best_split->out.count >= min_count)
                results.future_children.push_back(control_plane.enqueue([this, data, out_set = results.best_split->out.set, pop_name]() {
                    return do_Pursuit(data, out_set, pop_name); }));
        }
    }    

    return results;
}

void Pursuit_Results::wait() noexcept
{
    for (auto &child : future_children)
        children.push_back(child.get());
    for (auto &child : children)
        child.wait();
}    