// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// hello_world
void hello_world();
RcppExport SEXP _wdnet_hello_world() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    hello_world();
    return R_NilValue;
END_RCPP
}
// fx
double fx(arma::colvec x, arma::mat Y, arma::colvec z);
RcppExport SEXP _wdnet_fx(SEXP xSEXP, SEXP YSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(fx(x, Y, z));
    return rcpp_result_gen;
END_RCPP
}
// directed_rewire_cpp
Rcpp::List directed_rewire_cpp(int iteration, int nattempts, arma::uvec targetNode, arma::vec sourceOut, arma::vec sourceIn, arma::vec targetOut, arma::vec targetIn, arma::vec r_sourceOut, arma::vec r_sourceIn, arma::vec r_targetOut, arma::vec r_targetIn, arma::uvec index_s, arma::uvec index_t, arma::mat eta, bool rewireTrack);
RcppExport SEXP _wdnet_directed_rewire_cpp(SEXP iterationSEXP, SEXP nattemptsSEXP, SEXP targetNodeSEXP, SEXP sourceOutSEXP, SEXP sourceInSEXP, SEXP targetOutSEXP, SEXP targetInSEXP, SEXP r_sourceOutSEXP, SEXP r_sourceInSEXP, SEXP r_targetOutSEXP, SEXP r_targetInSEXP, SEXP index_sSEXP, SEXP index_tSEXP, SEXP etaSEXP, SEXP rewireTrackSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iteration(iterationSEXP);
    Rcpp::traits::input_parameter< int >::type nattempts(nattemptsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type targetNode(targetNodeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sourceOut(sourceOutSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sourceIn(sourceInSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type targetOut(targetOutSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type targetIn(targetInSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type r_sourceOut(r_sourceOutSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type r_sourceIn(r_sourceInSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type r_targetOut(r_targetOutSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type r_targetIn(r_targetInSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type index_s(index_sSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type index_t(index_tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< bool >::type rewireTrack(rewireTrackSEXP);
    rcpp_result_gen = Rcpp::wrap(directed_rewire_cpp(iteration, nattempts, targetNode, sourceOut, sourceIn, targetOut, targetIn, r_sourceOut, r_sourceIn, r_targetOut, r_targetIn, index_s, index_t, eta, rewireTrack));
    return rcpp_result_gen;
END_RCPP
}
// undirected_rewire_cpp
Rcpp::List undirected_rewire_cpp(int iteration, int nattempts, Rcpp::IntegerVector node1, Rcpp::IntegerVector node2, arma::vec degree1, arma::vec degree2, arma::vec index1, arma::vec index2, arma::mat e, bool rewireTrack);
RcppExport SEXP _wdnet_undirected_rewire_cpp(SEXP iterationSEXP, SEXP nattemptsSEXP, SEXP node1SEXP, SEXP node2SEXP, SEXP degree1SEXP, SEXP degree2SEXP, SEXP index1SEXP, SEXP index2SEXP, SEXP eSEXP, SEXP rewireTrackSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iteration(iterationSEXP);
    Rcpp::traits::input_parameter< int >::type nattempts(nattemptsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type node1(node1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type node2(node2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type degree1(degree1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type degree2(degree2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type index1(index1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type index2(index2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type e(eSEXP);
    Rcpp::traits::input_parameter< bool >::type rewireTrack(rewireTrackSEXP);
    rcpp_result_gen = Rcpp::wrap(undirected_rewire_cpp(iteration, nattempts, node1, node2, degree1, degree2, index1, index2, e, rewireTrack));
    return rcpp_result_gen;
END_RCPP
}
// findNode_cpp
arma::vec findNode_cpp(arma::vec nodes, arma::vec edges, arma::vec index);
RcppExport SEXP _wdnet_findNode_cpp(SEXP nodesSEXP, SEXP edgesSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type edges(edgesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(findNode_cpp(nodes, edges, index));
    return rcpp_result_gen;
END_RCPP
}
// nodeStrength_cpp
Rcpp::List nodeStrength_cpp(arma::vec startNode, arma::vec endNode, arma::vec weight, int nNodes, bool weighted);
RcppExport SEXP _wdnet_nodeStrength_cpp(SEXP startNodeSEXP, SEXP endNodeSEXP, SEXP weightSEXP, SEXP nNodesSEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type startNode(startNodeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type endNode(endNodeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type nNodes(nNodesSEXP);
    Rcpp::traits::input_parameter< bool >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(nodeStrength_cpp(startNode, endNode, weight, nNodes, weighted));
    return rcpp_result_gen;
END_RCPP
}
// sampleNode_cpp
arma::vec sampleNode_cpp(arma::vec totalNode);
RcppExport SEXP _wdnet_sampleNode_cpp(SEXP totalNodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type totalNode(totalNodeSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleNode_cpp(totalNode));
    return rcpp_result_gen;
END_RCPP
}
// rpanet_cpp
Rcpp::List rpanet_cpp(arma::vec startNode, arma::vec endNode, arma::vec scenario, int nNodes, int nEdges, double delta_out, double delta_in);
RcppExport SEXP _wdnet_rpanet_cpp(SEXP startNodeSEXP, SEXP endNodeSEXP, SEXP scenarioSEXP, SEXP nNodesSEXP, SEXP nEdgesSEXP, SEXP delta_outSEXP, SEXP delta_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type startNode(startNodeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type endNode(endNodeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type scenario(scenarioSEXP);
    Rcpp::traits::input_parameter< int >::type nNodes(nNodesSEXP);
    Rcpp::traits::input_parameter< int >::type nEdges(nEdgesSEXP);
    Rcpp::traits::input_parameter< double >::type delta_out(delta_outSEXP);
    Rcpp::traits::input_parameter< double >::type delta_in(delta_inSEXP);
    rcpp_result_gen = Rcpp::wrap(rpanet_cpp(startNode, endNode, scenario, nNodes, nEdges, delta_out, delta_in));
    return rcpp_result_gen;
END_RCPP
}
// sampleNode_simple_cpp
int sampleNode_simple_cpp(int tnode, double sumstrength, arma::vec strength, double delta);
RcppExport SEXP _wdnet_sampleNode_simple_cpp(SEXP tnodeSEXP, SEXP sumstrengthSEXP, SEXP strengthSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type tnode(tnodeSEXP);
    Rcpp::traits::input_parameter< double >::type sumstrength(sumstrengthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type strength(strengthSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleNode_simple_cpp(tnode, sumstrength, strength, delta));
    return rcpp_result_gen;
END_RCPP
}
// rpanet_simple_cpp
Rcpp::List rpanet_simple_cpp(int nsteps, arma::vec control, arma::vec m, arma::vec w, arma::vec outstrength, arma::vec instrength, double sumstrength, int nnode);
RcppExport SEXP _wdnet_rpanet_simple_cpp(SEXP nstepsSEXP, SEXP controlSEXP, SEXP mSEXP, SEXP wSEXP, SEXP outstrengthSEXP, SEXP instrengthSEXP, SEXP sumstrengthSEXP, SEXP nnodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type control(controlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type outstrength(outstrengthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type instrength(instrengthSEXP);
    Rcpp::traits::input_parameter< double >::type sumstrength(sumstrengthSEXP);
    Rcpp::traits::input_parameter< int >::type nnode(nnodeSEXP);
    rcpp_result_gen = Rcpp::wrap(rpanet_simple_cpp(nsteps, control, m, w, outstrength, instrength, sumstrength, nnode));
    return rcpp_result_gen;
END_RCPP
}
