
#' initConfiguration initializes E2Predictor configuration
#'
#' @param netMHCpath path to netMHC 4.0 (only necessary if new allele predictions are needed)
#' @param working.path by default the current working path
#' @param hla_alleles list of cell.line - alleles
#'
#' @export
#'
initConfiguration <- function(
    #path to netMHC 4.0, if netMHC 4.0 is available
    netMHCpath = "/Users/napedro/external_sources/netMHC-4.0/netMHC", #NULL,
    working.path = getwd(),
    hla_alleles = list( JY =          c("A*02:01","B*07:02","C*07:02"),
                      LCLC.103H =   c("A*03:01","B*07:02","C*07:02"),
                      COLO.699.N =  c("A*02:01","B*07:02","C*07:02","C*06:02"),
                      MEL.526 =     c("A*02:01","B*15:01","C*07:02"),
                      HEK293 =      c("A*03:01","B*07:02","C*07:02"),
                      MEL.HO =      c("A*03:01","B*07:02","B*14:02","C*07:02","C*08:05"),
                      NIH.OVCAR.3 = c("A*02:01","A*29:02","B*07:02","B*58:01","C*07:02","C*07:19"),
                      SK.MEL.37 =   c("A*02:01","B*07:02","B*56:01","C*07:02","C*03:04"),
                      K562A2 =      c("A*02:01","C*05:01","C*03:04")
    )
){
    # list argument names
    argNames = ls()
    # collect arguments and their values in a list
    E2Predictor.Config <<- as.list( sapply( argNames, function(n) get(n) ) )
    # process configuration
    #E2Predictor.setDataRootFolder(E2Predictor.Config$DataRootFolder, createSubfolders=F)
    ################################################################################
    #E2Predictor.processConfig()

}


#' changeConfigutation allows changes at E2Predictor configuration
#'
#' @param netMHCpath path to netMHC 4.0 (only necessary if new allele predictions are needed)
#' @param working.path by default the current working path
#' @param hla_alleles list of cell.line - alleles
#'
#' @export
#'
changeConfiguration <- function(
    netMHCpath = E2Predictor.Config$netMHCpath,
    working.path = E2Predictor.Config$working.path,
    hla_alleles = E2Predictor.Config$hla_alleles
){
    # list argument names
    argNames = ls()
    # collect arguments and their values in a list
    E2Predictor.Config <<- as.list(sapply( argNames, function(n) get(n) ))
    # process configuration
    #E2Predictor.setDataRootFolder(E2Predictor.Config$DataRootFolder, createSubfolders=F)
    ################################################################################
    #E2Predictor.processConfig()
}
