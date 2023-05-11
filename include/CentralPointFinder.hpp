
#include "Initializer.hpp"


class CentralPointFinder : public Initializer {
    public:

    CentralPointFinder(const double ss, const double ts, const double te, 
                       const double err_term, const double grad_lim) : 
                       ss_(ss), ts_(ts), te_(te), err_term_(err_term), grad_lim_(grad_lim), Initializer() {}

    virtual VectorXd getInitialPoint(MatrixXd A, VectorXd b);


    protected:
        const double ss_;
        const double ts_;
        const double te_;
        const double err_term_;
        const double grad_lim_;

};