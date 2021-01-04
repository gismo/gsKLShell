/** @file gsMaterialMatrix.h

    @brief Provides material matrices for the thin shell class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

namespace gismo
{

template<class T> class gsMaterialMatrixBase;


/** @brief Assembles system matrices for thin shell linear and nonlinear elasticity problems.

    \tparam T coefficient type

    \ingroup gsThinShell
*/

/*
Desired template parameters
- compressible                      bool    COM
- material model                    int     MAT
- output type                       int     OUT / VEC
- thickness integration method      int     INT
- geometric dimension               int     DIM
*/

template <short_t dim, class T, short_t mat, bool comp>
class gsMaterialMatrix : public gsMaterialMatrixBase<T>
{
public:
    /** @brief Constructor of the assembler object.

        \param[in] ...
        \param[in] ...

    */
    /// Default empty constructor
    gsMaterialMatrix() { }

    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars);

    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars);

    gsMaterialMatrix(   const gsFunctionSet<T> & mp,
                        const gsFunctionSet<T> & mp_def,
                        const gsFunction<T> & thickness,
                        const std::vector<gsFunction<T> *> &pars,
                        const gsFunction<T> & Density);

    // Laminates without deformed
    gsMaterialMatrix(   const gsFunctionSet<T>            & mp,
                        const std::vector<T>                thickness,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T>              & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T>                phi);
    gsMaterialMatrix(   const gsFunctionSet<T>            & mp,
                        const std::vector<T>                thickness,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T>              & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T>                phi,
                        const std::vector<T>                density);

    // Laminates with deformed
    gsMaterialMatrix(   const gsFunctionSet<T>            & mp,
                        const gsFunctionSet<T>            & mp_def,
                        const std::vector<T>                thickness,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T>              & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T>                phi);
    gsMaterialMatrix(   const gsFunctionSet<T>            & mp,
                        const gsFunctionSet<T>            & mp_def,
                        const std::vector<T>                thickness,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T>              & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T>                phi,
                        const std::vector<T>                density);


    /// @brief Returns the list of default options for assembly
    gsOptionList & options() {return m_options;}
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt



    /// Shared pointer for gsMaterialMatrix
    typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

    /// Unique pointer for gsMaterialMatrix
    typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // CORRECT???
    // ~gsMaterialMatrix() { delete m_piece; }

    short_t domainDim() const;

    // template OUT
    short_t targetDim() const;

    const gsFunction<T> & piece(const index_t p) const
    {
        // delete m_piece;
        // m_piece = new gsMaterialMatrix(m_patches->piece(k), *m_thickness, *m_YoungsModulus, *m_PoissonRatio);
        m_piece = new gsMaterialMatrix(*this);
        m_piece->setPatch(p);
        return *m_piece;
    }

    ~gsMaterialMatrix() { delete m_piece; }

    void setPatch(index_t p) {m_pIndex = p; }

    // template OUT INT
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // template COMPOSITE
    void eval_into_dens(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // template COM
    void eval_into_stretch(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // template COM
    void eval_into_stretchdir(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    // template OUT VEC INT COMPOSITE
    void eval_into_NP(const gsMatrix<T>& u, gsMatrix<T>& result) const;

protected:
    template<bool _com> inline
    typename std::enable_if<_com, void>::type eval_into_stretch_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    template<bool _com> inline
    typename std::enable_if<!_com, void>::type eval_into_stretch_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    template<bool _com> inline
    typename std::enable_if<_com, void>::type eval_into_stretchdir_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    template<bool _com> inline
    typename std::enable_if<!_com, void>::type eval_into_stretchdir_impl(const gsMatrix<T>& u, gsMatrix<T>& result) const;

public:
    // template<short_t num=0>
    // void makeMatrix()                   { m_outputType=2; m_output = num;}

    template<short_t num=0>
    gsMaterialMatrix * makeMatrix()
    {
        gsMaterialMatrix * tmp = new gsMaterialMatrix(*this);
        tmp->m_outputType = 2;
        tmp->m_output = num;
        return tmp;
    }


    template<short_t num=0>
    gsMaterialMatrix * makeVector()
    {
        gsMaterialMatrix * tmp = new gsMaterialMatrix(*this);
        tmp->m_outputType = 1;
        tmp->m_output = num;
        return tmp;
    }

    template<short_t num=0>
    gsMaterialMatrix * makePrincipleStress()
    {
        gsMaterialMatrix * tmp = new gsMaterialMatrix(*this);
        tmp->m_outputType = 10;
        tmp->m_output = num;
        return tmp;
    }

    gsMaterialMatrix * makeDensity()
    {
        gsMaterialMatrix * tmp = new gsMaterialMatrix(*this);
        tmp->m_outputType = 0;
        return tmp;
    }

    gsMaterialMatrix * makeStretch()
    {
        gsMaterialMatrix * tmp = new gsMaterialMatrix(*this);
        tmp->m_outputType = 9;
        return tmp;
    }

    gsMaterialMatrix * makeDirections()
    {
        gsMaterialMatrix * tmp = new gsMaterialMatrix(*this);
        tmp->m_outputType = 11;
        return tmp;
    }

    void setParameters(const std::vector<gsFunction<T>*> &pars)
    {
        m_pars = pars;
        m_numPars = m_pars.size();
    }

    void info() const;


protected:
    void initialize();
    void defaultOptions();
    void getOptions() const;

    // template COM MAT
    gsMatrix<T> eval3D(const index_t i, const gsMatrix<T>& z) const;
    gsMatrix<T> eval3D(const index_t i) const;
    // template OUT VEC
    gsMatrix<T> eval_Compressible(const index_t i, const gsMatrix<T>& z) const;
    // template OUT VEC
    gsMatrix<T> eval_Incompressible(const index_t i, const gsMatrix<T>& z) const;

    private:

        template<short_t _mat, bool _com> inline
        typename std::enable_if<_mat==0, gsMatrix<T>>::type eval3D_impl(const index_t i, const gsMatrix<T>& z) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<_com && !(_mat==0), gsMatrix<T>>::type eval3D_impl(const index_t i, const gsMatrix<T>& z) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<!_com && !(_mat==0), gsMatrix<T>>::type eval3D_impl(const index_t i, const gsMatrix<T>& z) const;

    protected:



    // template MAT
    T Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l) const;
    // template MAT
    T Cijkl3D(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    // template MAT
    T Cijkl  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    // template MAT
    T Sij    (const index_t i, const index_t j) const;
    // template MAT
    T Sij    (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    T Sii    (const index_t i) const;
    T Sii    (const index_t i, const gsMatrix<T> & c) const;

    private:

        // Incompressible
        // template<short_t _mat> inline
        // typename std::enable_if<_mat==0, T>::type Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==0, T>::type Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==2, T>::type Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==3, T>::type Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==5, T>::type Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
        {GISMO_NO_IMPLEMENTATION};
        template<short_t _mat> inline
        typename std::enable_if<(_mat >= 10) && (_mat < 20) , T>::type Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
        template<short_t _mat> inline
        typename std::enable_if<(_mat >= 20) && (_mat < 30) , T>::type Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;

        // Compressible
        template<short_t _mat> inline
        typename std::enable_if<(_mat >= 10) && (_mat < 20),   T>::type Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<!((_mat >= 10) && (_mat < 20)),T>::type Cijkl_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

        template<short_t _mat> inline
        typename std::enable_if<_mat==0, T>::type Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==2, T>::type Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==3, T>::type Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==5, T>::type Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<(_mat >= 10) && (_mat < 20) , T>::type Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<(_mat >= 20) && (_mat < 30) , T>::type Cijkl3D_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

        template<short_t _mat> inline
        typename std::enable_if<_mat==0, T>::type Sij_impl(const index_t i, const index_t j) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==2, T>::type Sij_impl(const index_t i, const index_t j) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==3, T>::type Sij_impl(const index_t i, const index_t j) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==5, T>::type Sij_impl(const index_t i, const index_t j) const
        {GISMO_NO_IMPLEMENTATION};
        template<short_t _mat> inline
        typename std::enable_if<(_mat >= 10) && (_mat < 20) , T>::type Sij_impl(const index_t i, const index_t j) const;
        template<short_t _mat> inline
        typename std::enable_if<(_mat >= 20) && (_mat < 30) , T>::type Sij_impl(const index_t i, const index_t j) const;

        template<short_t _mat> inline
        typename std::enable_if<_mat==0, T>::type Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
        {GISMO_NO_IMPLEMENTATION};
        template<short_t _mat> inline
        typename std::enable_if<_mat==2, T>::type Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==3, T>::type Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==5, T>::type Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<(_mat >= 10) && (_mat < 20) , T>::type Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<(_mat >= 20) && (_mat < 30) , T>::type Sij_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    protected:

    // template MAT
    T dPsi   (const index_t i, const index_t j) const;
    // template MAT
    T dPsi   (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
    // template MAT
    T d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l) const;
    // template MAT
    T d2Psi  (const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    T dI_1   (const index_t i, const index_t j) const;
    T dI_2   (const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;

    private:
        /*
            Note: when adding a new implementation, make sure to add the condition in the 'other' implementation : !(_mat==xx || ...)
        */
        template<short_t _mat> inline
        typename std::enable_if<_mat==22, T>::type dPsi_impl(const index_t i, const index_t j) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==23, T>::type dPsi_impl(const index_t i, const index_t j) const;
        // other
        template<short_t _mat> inline
        typename std::enable_if<!(_mat==22 || _mat==23), T>::type dPsi_impl(const index_t i, const index_t j) const
        {GISMO_NO_IMPLEMENTATION};
        // add other

        template<short_t _mat> inline
        typename std::enable_if<_mat==22, T>::type dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==23, T>::type dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==25, T>::type dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        // other
        template<short_t _mat> inline
        typename std::enable_if<!(_mat==22 || _mat==23 || _mat==25), T>::type dPsi_impl(const index_t i, const index_t j, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
        {GISMO_NO_IMPLEMENTATION};

        template<short_t _mat> inline
        typename std::enable_if<_mat==22, T>::type d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==23, T>::type d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const;
        // other
        template<short_t _mat> inline
        typename std::enable_if<!(_mat==22 || _mat==23), T>::type d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l) const
        {GISMO_NO_IMPLEMENTATION};

        template<short_t _mat> inline
        typename std::enable_if<_mat==22, T>::type d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==23, T>::type d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat==25, T>::type d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const;
        // other
        template<short_t _mat> inline
        typename std::enable_if<!(_mat==22 || _mat==23 || _mat==25), T>::type d2Psi_impl(const index_t i, const index_t j, const index_t k, const index_t l, const gsMatrix<T> & c, const gsMatrix<T> & cinv) const
        {GISMO_NO_IMPLEMENTATION};
    protected:

    // Stretch based formulation
    // template MAT COM
    T dPsi_da   (const index_t a) const;
    // template MAT COM
    T d2Psi_dab (const index_t a, const index_t b) const;
    T dJ_da     (const index_t a) const;
    T d2J_dab   (const index_t a, const index_t b) const;
    T p()                                          const;
    T dp_da     (const index_t a) const;
    // template COM
    T Sa        (const index_t a) const;
    // template COM
    T dSa_db    (const index_t a, const index_t b) const;
    // template COM
    T Cabcd     (const index_t a, const index_t b, const index_t c, const index_t d) const;

    private:
        // ----------------------------------------------------------------------------------

        template<short_t _mat, bool _com> inline
        typename std::enable_if<_com && (_mat==2 || _mat==12 || _mat==22), T>::type dPsi_da_impl(const index_t a) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<!_com && (_mat==2 || _mat==12 || _mat==22), T>::type dPsi_da_impl(const index_t a) const;

        template<short_t _mat, bool _com> inline
        typename std::enable_if<_com && (_mat==3 || _mat==13 || _mat==23), T>::type dPsi_da_impl(const index_t a) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<!_com && (_mat==3 || _mat==13 || _mat==23), T>::type dPsi_da_impl(const index_t a) const;

        template<short_t _mat, bool _com> inline
        typename std::enable_if<_com && _mat==14, T>::type dPsi_da_impl(const index_t a) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<!_com && _mat==14, T>::type dPsi_da_impl(const index_t a) const;

        template<short_t _mat, bool _com> inline
        typename std::enable_if<_com && (_mat==5 || _mat==15 || _mat==25), T>::type dPsi_da_impl(const index_t a) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<!_com && (_mat==5 || _mat==15 || _mat==25), T>::type dPsi_da_impl(const index_t a) const
        {GISMO_NO_IMPLEMENTATION};

        // other
        template<short_t _mat, bool _com> inline
        typename std::enable_if<
                                !(_mat==2 || _mat==12 || _mat==22) &&
                                !(_mat==3 || _mat==13 || _mat==23) &&
                                !(_mat==14) &&
                                !(_mat==5 || _mat==15 || _mat==25)
                                                                    , T>::type dPsi_da_impl(const index_t a) const
        {GISMO_NO_IMPLEMENTATION}

        // ----------------------------------------------------------------------------------

        template<short_t _mat, bool _com> inline
        typename std::enable_if<_com && (_mat==2 || _mat==12 || _mat==22), T>::type d2Psi_dab_impl(const index_t a, const index_t b) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<!_com && (_mat==2 || _mat==12 || _mat==22), T>::type d2Psi_dab_impl(const index_t a, const index_t b) const;

        template<short_t _mat, bool _com> inline
        typename std::enable_if<_com && (_mat==3 || _mat==13 || _mat==23), T>::type d2Psi_dab_impl(const index_t a, const index_t b) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<!_com && (_mat==3 || _mat==13 || _mat==23), T>::type d2Psi_dab_impl(const index_t a, const index_t b) const;

        template<short_t _mat, bool _com> inline
        typename std::enable_if<_com && _mat==14, T>::type d2Psi_dab_impl(const index_t a, const index_t b) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<!_com && _mat==14, T>::type d2Psi_dab_impl(const index_t a, const index_t b) const;

        template<short_t _mat, bool _com> inline
        typename std::enable_if<_com && (_mat==5 || _mat==15 || _mat==25), T>::type d2Psi_dab_impl(const index_t a, const index_t b) const;
        template<short_t _mat, bool _com> inline
        typename std::enable_if<!_com && (_mat==5 || _mat==15 || _mat==25), T>::type d2Psi_dab_impl(const index_t a, const index_t b) const
        {GISMO_NO_IMPLEMENTATION};

        // other
        template<short_t _mat, bool _com> inline
        typename std::enable_if<
                                !(_mat==2 || _mat==12 || _mat==22) &&
                                !(_mat==3 || _mat==13 || _mat==23) &&
                                !(_mat==14) &&
                                !(_mat==5 || _mat==15 || _mat==25)
                                                                    , T>::type d2Psi_dab_impl(const index_t a, const index_t b) const
        {GISMO_NO_IMPLEMENTATION}

        // ----------------------------------------------------------------------------------

        template<bool _com> inline
        typename std::enable_if<_com , T>::type Sa_impl(const index_t a) const;
        template<bool _com> inline
        typename std::enable_if<!_com, T>::type Sa_impl(const index_t a) const;

        // ----------------------------------------------------------------------------------

        template<bool _com> inline
        typename std::enable_if<_com , T>::type dSa_db_impl(const index_t a, const index_t b) const;
        template<bool _com> inline
        typename std::enable_if<!_com, T>::type dSa_db_impl(const index_t a, const index_t b) const;

        // ----------------------------------------------------------------------------------

        template<bool _com> inline
        typename std::enable_if<_com , T>::type Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const;
        template<bool _com> inline
        typename std::enable_if<!_com, T>::type Cabcd_impl(const index_t a, const index_t b, const index_t c, const index_t d) const;

        // ----------------------------------------------------------------------------------

    protected:

    // T Sij  (const index_t i, const index_t j, gsMatrix<T> & cinv) const;

    // void eval_Linear(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // void eval_Incompressible(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    // void eval_Compressible(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    gsMatrix<T> integrateZ(const gsMatrix<T>& u) const;
    gsMatrix<T> multiplyZ (const gsMatrix<T>& u) const;

    // template MAT
    void computePoints(const gsMatrix<T> & u, bool deformed=true) const;
    // template DIM
    void computeMetricDeformed() const;
    // template DIM
    void computeMetricUndeformed() const;
    // template DIM
    void getMetric(index_t k, T z) const;
    // template DIM
    void getMetricDeformed(index_t k, T z) const;
    // template DIM
    void getMetricUndeformed(index_t k, T z) const;
    // template DIM
    void computeStretch(const gsMatrix<T> & C ) const;

    // template DIM
    std::pair<gsVector<T>,gsMatrix<T>> evalStretch(const gsMatrix<T> & C ) const;

    private:
        template<short_t _mat> inline
        typename std::enable_if<_mat==14, void>::type computePoints_impl(const gsMatrix<T> & u, bool deformed) const;
        template<short_t _mat> inline
        typename std::enable_if<_mat!=14, void>::type computePoints_impl(const gsMatrix<T> & u, bool deformed) const;

        template<short_t _dim> inline
        typename std::enable_if<_dim==2, void>::type computeMetricDeformed_impl() const;
        template<short_t _dim> inline
        typename std::enable_if<_dim==3, void>::type computeMetricDeformed_impl() const;

        template<short_t _dim> inline
        typename std::enable_if<_dim==2, void>::type computeMetricUndeformed_impl() const;
        template<short_t _dim> inline
        typename std::enable_if<_dim==3, void>::type computeMetricUndeformed_impl() const;

        template<short_t _dim> inline
        typename std::enable_if<_dim==2, void>::type getMetric_impl(index_t k, T z) const;
        template<short_t _dim> inline
        typename std::enable_if<_dim==3, void>::type getMetric_impl(index_t k, T z) const;

        template<short_t _dim> inline
        typename std::enable_if<_dim==2, void>::type getMetricDeformed_impl(index_t k, T z) const;
        template<short_t _dim> inline
        typename std::enable_if<_dim==3, void>::type getMetricDeformed_impl(index_t k, T z) const;

        template<short_t _dim> inline
        typename std::enable_if<_dim==2, void>::type getMetricUndeformed_impl(index_t k, T z) const;
        template<short_t _dim> inline
        typename std::enable_if<_dim==3, void>::type getMetricUndeformed_impl(index_t k, T z) const;

        template<short_t _dim> inline
        typename std::enable_if<_dim==2, void>::type computeStretch_impl(const gsMatrix<T> & C ) const;
        template<short_t _dim> inline
        typename std::enable_if<_dim==3, void>::type computeStretch_impl(const gsMatrix<T> & C ) const;

        template<short_t _dim> inline
        typename std::enable_if<_dim==2, std::pair<gsVector<T>,gsMatrix<T>>>::type evalStretch_impl(const gsMatrix<T> & C ) const;
        template<short_t _dim> inline
        typename std::enable_if<_dim==3, std::pair<gsVector<T>,gsMatrix<T>>>::type evalStretch_impl(const gsMatrix<T> & C ) const;


protected:
    // general
    index_t m_pIndex;
    index_t m_numPars; // how many parameters for the material model?
    mutable int m_moment;
    mutable gsMatrix<T> m_result;

    // constructor
    const gsFunctionSet<T> * m_patches;
    const gsFunctionSet<T> * m_defpatches;
    const gsFunction<T> * m_thickness;
    mutable std::vector<gsFunction<T>* > m_pars;
    const gsFunction<T> * m_density;

    // Linear material matrix
    mutable gsMapData<T> m_map, m_map_def;
    mutable gsMatrix<T,2,2> m_Acov_ori, m_Acon_ori, m_Acov_def, m_Acon_def, m_Bcov_ori, m_Bcon_ori, m_Bcov_def, m_Bcon_def;
    mutable gsMatrix<T,3,2> m_acov_ori, m_acon_ori, m_acov_def, m_acon_def, m_ncov_ori, m_ncov_def;
    mutable gsMatrix<T> m_Acov_ori_mat, m_Acon_ori_mat, m_Acov_def_mat, m_Acon_def_mat, m_Bcov_ori_mat, m_Bcov_def_mat;
    mutable gsMatrix<T> m_acov_ori_mat, m_acon_ori_mat, m_acov_def_mat, m_acon_def_mat, m_ncov_ori_mat, m_ncov_def_mat;
    mutable gsMatrix<T> m_Emat,m_Nmat,m_Tmat,m_rhomat;
    mutable real_t m_lambda, m_mu, m_Cconstant;


    // Composite material matrix
    const std::vector<std::pair<T,T>>   m_YoungsModuli;
    const std::vector<T>                m_ShearModuli;
    const std::vector<std::pair<T,T>>   m_PoissonRatios;
    const std::vector<T>                m_thickValues;
    const std::vector<T>                m_phis;
    const std::vector<T>                m_densities;
    mutable gsVector<T>                 m_e1, m_e2, m_ac1, m_ac2, m_a1, m_a2;
    mutable gsMatrix<T>                 m_covBasis, m_covMetric, m_conBasis, m_conMetric;
    mutable gsMatrix<T>                 m_stretches, m_stretchvec;
    mutable T                           m_E1, m_E2, m_G12, m_nu12, m_nu21, m_t, m_t_tot,
                                        m_t_temp, m_z, m_z_mid, m_phi, m_rho;
    mutable gsMatrix<T>                 m_Dmat, m_Transform;

    // Compressible material matrix
    mutable gsMatrix<T>                 m_deriv2_def, m_deriv2_ori;
    mutable gsMatrix<T,3,3>             m_Gcov_ori, m_Gcon_ori, m_Gcov_def, m_Gcon_def, m_Gcov_ori_L, m_Gcov_def_L;
    mutable gsMatrix<T,3,3>             m_gcov_ori, m_gcov_def,m_gcon_ori, m_gcon_def;
    mutable gsMatrix<T>                 m_parmat;
    mutable gsVector<T>                 m_parvals;
    mutable T                           m_J0, m_J0_sq, m_J, m_J_sq, m_Tval;
    // integrateZ
    mutable gsMatrix<T> m_points2D, m_points3D, m_evalPoints;
    // mutable gsMatrix<T> m_quNodes;
    // mutable gsVector<T> m_quWeights;
    mutable gsGaussRule<T> m_gauss;
    mutable index_t m_numGauss;
    mutable T m_tHalf;


    mutable gsOptionList m_options;

    mutable gsMaterialMatrix<dim,T,mat,comp> * m_piece; // todo: improve the way pieces are accessed

    /// @brief Specifies the material law to use
    struct material_law
    {
        enum type
        {
            SvK_Isotropic = 0,          /// Psi = ........ S = 2*mu*E + lambda*tr(E)*I
            SvK_Orthotropic = 1,        /// Psi = ........ S = 2*mu*E + lambda*tr(E)*I
            NHK = 2,       /// Psi = ........ S = lambda*ln(J)*C^-1 + mu*(I-C^-1)
            MR = 3        /// Psi = ........ S = lambda*ln(J)*C^-1 + mu*(I-C^-1)
        };
    };
    /// @brief Specifies (in)compressibility
    struct compressibility
    {
        enum type
        {
            incompressible = 0,
            compressible = 1
        };
    };
    /// @brief Specifies (in)compressibility
    struct integration
    {
        enum type
        {
            DD = 0,
            AP = 1,
            NP = 2
        };
    };

    mutable index_t m_material;
    mutable bool m_compressible;
    mutable int m_compFun;
    mutable int m_outputType, m_output;
    mutable int m_integration;



};

template <class T>
class gsMaterialMatrixBase : public gsFunction<T>
{
public:
    virtual ~gsMaterialMatrixBase() {};

    // GISMO_CLONE_FUNCTION(gsMaterialMatrixBase)


    inline virtual gsOptionList & options() = 0;
    inline virtual void setOptions(gsOptionList opt) = 0;

    inline virtual short_t domainDim() const = 0;

    // template OUT
    inline virtual short_t targetDim() const = 0;

    inline virtual const gsFunction<T> & piece(const index_t p) const = 0;

    inline virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;
    inline virtual void eval_into_dens(const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;
    inline virtual void eval_into_stretch(const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;
    inline virtual void eval_into_stretchdir(const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;
    inline virtual void eval_into_NP(const gsMatrix<T>& u, gsMatrix<T>& result) const = 0;

    // template<short_t num=0>
    inline virtual gsMaterialMatrixBase * makeMatrix() = 0;
    // template<short_t num=0>
    inline virtual gsMaterialMatrixBase * makeVector() = 0;
    // template<short_t num=0>
    inline virtual gsMaterialMatrixBase * makePrincipleStress() = 0;

    inline virtual gsMaterialMatrixBase * makeDensity() = 0;
    inline virtual gsMaterialMatrixBase * makeStretch() = 0;
    inline virtual gsMaterialMatrixBase * makeDirections() = 0;

    inline virtual void setParameters(const std::vector<gsFunction<T>*> &pars) =0;
    inline virtual void info() const = 0;
};

} // namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMaterialMatrix.hpp)
#endif


// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
// template <class T>
// class gsMaterialMatrix : public gismo::gsFunction<T>
// {
//   // Computes the material matrix for different material models
//   //
// protected:
//     const gsFunctionSet<T> * _mp;
//     const gsFunction<T> * _YoungsModulus;
//     const gsFunction<T> * _PoissonRatio;
//     mutable gsMapData<T> _tmp;
//     mutable gsMatrix<real_t,3,3> F0;
//     mutable gsMatrix<T> Emat,Nmat;
//     mutable real_t lambda, mu, E, nu, C_constant;

// public:
//     /// Shared pointer for gsMaterialMatrix
//     typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

//     /// Unique pointer for gsMaterialMatrix
//     typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW

//     gsMaterialMatrix(const gsFunctionSet<T> & mp, const gsFunction<T> & YoungsModulus,
//                    const gsFunction<T> & PoissonRatio) :
//     _mp(&mp), _YoungsModulus(&YoungsModulus), _PoissonRatio(&PoissonRatio), _mm_piece(nullptr)
//     {
//         _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
//     }

//     ~gsMaterialMatrix() { delete _mm_piece; }

//     GISMO_CLONE_FUNCTION(gsMaterialMatrix)

//     short_t domainDim() const {return 2;}

//     short_t targetDim() const {return 9;}

//     mutable gsMaterialMatrix<T> * _mm_piece; // todo: improve the way pieces are accessed

//     const gsFunction<T> & piece(const index_t k) const
//     {
//         delete _mm_piece;
//         _mm_piece = new gsMaterialMatrix(_mp->piece(k), *_YoungsModulus, *_PoissonRatio);
//         return *_mm_piece;
//     }

//     //class .. matMatrix_z
//     // should contain eval_into(thickness variable)

//     // Input is parametric coordinates of the surface \a mp
//     void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
//     {
//         // NOTE 1: if the input \a u is considered to be in physical coordinates
//         // then we first need to invert the points to parameter space
//         // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
//         // otherwise we just use the input paramteric points
//         _tmp.points = u;

//         static_cast<const gsFunction<T>&>(_mp->piece(0)).computeMap(_tmp); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

//         // NOTE 2: in the case that parametric value is needed it suffices
//         // to evaluate Youngs modulus and Poisson's ratio at
//         // \a u instead of _tmp.values[0].
//         _YoungsModulus->eval_into(_tmp.values[0], Emat);
//         _PoissonRatio->eval_into(_tmp.values[0], Nmat);

//         result.resize( targetDim() , u.cols() );
//         for( index_t i=0; i< u.cols(); ++i )
//         {
//             gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);

//             F0.leftCols(2) = _tmp.jacobian(i);
//             F0.col(2)      = _tmp.normal(i).normalized();
//             F0 = F0.inverse();
//             F0 = F0 * F0.transpose(); //3x3

//             // Evaluate material properties on the quadrature point
//             E = Emat(0,i);
//             nu = Nmat(0,i);
//             lambda = E * nu / ( (1. + nu)*(1.-2.*nu)) ;
//             mu     = E / (2.*(1. + nu)) ;

//             C_constant = 2*lambda*mu/(lambda+2*mu);

//             C(0,0) = C_constant*F0(0,0)*F0(0,0) + 1*mu*(2*F0(0,0)*F0(0,0));
//             C(1,1) = C_constant*F0(1,1)*F0(1,1) + 1*mu*(2*F0(1,1)*F0(1,1));
//             C(2,2) = C_constant*F0(0,1)*F0(0,1) + 1*mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
//             C(1,0) =
//             C(0,1) = C_constant*F0(0,0)*F0(1,1) + 1*mu*(2*F0(0,1)*F0(0,1));
//             C(2,0) =
//             C(0,2) = C_constant*F0(0,0)*F0(0,1) + 1*mu*(2*F0(0,0)*F0(0,1));
//             C(2,1) = C(1,2) = C_constant*F0(0,1)*F0(1,1) + 1*mu*(2*F0(0,1)*F0(1,1));

//             //gsDebugVar(C);
//         }
//     }

//     // piece(k) --> for patch k

// }; //! [Include namespace]

