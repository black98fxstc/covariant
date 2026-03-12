template <unsigned Dimension, unsigned Grid>
class Covariant {
public:
    const unsigned short stride[Dimension];
    float weight[static_cast<size_t>(std::pow(Grid, Dimension))];
};