#include <cmath>
#include <vector>

template <unsigned Dimension, unsigned Grid>
class Covariant
{
public:
    static constexpr unsigned dimension = Dimension;
    static constexpr unsigned grid = Grid;
    size_t stride[Dimension];
    typedef std::vector<float> Event;

private:
    static constexpr unsigned points = Grid + 1;
    size_t _events = 0;
    std::vector<float> _weight;

public:
    bool event(Event& event) {
        size_t coordinate = 0;
        double rem[Dimension];
        for (unsigned i = 0; i < dimension; i++)
        {
            if (event[i] < 0.0f || event[i] >= 1.0f)
                return false;
            rem[i] = event[i] * grid;
            unsigned floor = static_cast<unsigned>(rem[i]);
            coordinate += floor * stride[i];
            rem[i] -= floor;
        }
        
        for (unsigned neighbor = 0; neighbor < 1 << dimension; neighbor++)
        {
            size_t offset = 0;
            double weight = 1.0;
            for (unsigned i = 0; i < Dimension; i++)
            {
                if (neighbor & (1 << i)) {
                    offset += stride[i];
                    weight *= rem[i];
                } else {
                    weight *= 1.0f - rem[i];
                }
            }
            _weight[coordinate + offset] += weight;
        }
        ++_events;
        return true;        
    }

    size_t events()
    {
        return _events;
    }

    Covariant()
        : _weight(static_cast<size_t>(std::pow(points, dimension)))
    {
        size_t s = 1;
        for (unsigned i = 0; i < dimension; i++)
        {
            stride[i] = s;
            s *= points;
        }
    }

    size_t event_count() const
    {
        return _events;
    }
};
