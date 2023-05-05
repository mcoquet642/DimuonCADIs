
#ifndef buildCharmoniaMassModel_C
#define buildCharmoniaMassModel_C

#include "Utilities/initClasses.h"
#include "Utilities/RooExtCBShape.h"

void setMassDefaultParameters(map<string, string> &parIni, double numEntries);
bool addSignalMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni); 
bool addBackgroundMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni);


bool buildCharmoniaMassModel(RooWorkspace& ws, struct CharmModel model, map<string, string>  parIni, 
                             bool doConstrFit,            // Do constrained fit
                             bool incBkg,                 // Include background model
                             bool incJpsi,                // Include Jpsi model
                             double  numEntries = 300000. // Number of entries in the dataset
                             )
{

  // If the initial parameters are empty, set defaul parameter values
  setMassDefaultParameters(parIni, numEntries);


  // C r e a t e   m o d e l  

  if (incJpsi) {
    if(!addSignalMassModel(ws, "Jpsi", model.Jpsi.Mass, parIni)) { cout << "[ERROR] Adding Jpsi Mass Model failed" << endl; return false; }
  }
  if (incBkg) {
    if(!addBackgroundMassModel(ws, "Bkg", model.Bkg.Mass, parIni)) { cout << "[ERROR] Adding Background Mass Model failed" << endl; return false; }
  }
  // Constraint PDFs
  if (doConstrFit) //FIXME: hardcoded values should be moved to input files
  {
      ws.factory(Form("Gaussian::sigmaAlphaConstr(%s,RooConstVar(%f),RooConstVar(%f))",Form("alpha_Jpsi_%s", "PP"), ws.var(Form("alpha_Jpsi_%s", "PP"))->getValV(), 0.16*ws.var(Form("alpha_Jpsi_%s", "PP"))->getValV()));
      ws.factory(Form("Gaussian::sigmaNConstr(%s,RooConstVar(%f),RooConstVar(%f))",Form("n_Jpsi_%s", "PP"), ws.var(Form("n_Jpsi_%s", "PP"))->getValV(), 0.21*ws.var(Form("n_Jpsi_%s", "PP"))->getValV()));
  }
  
  // Total PDF
  string pdfType = "pdfMASS";
  string pdfName = Form("%s_Tot_%s", pdfType.c_str(), "PP");

  RooArgList pdfList;
  if (incJpsi) { pdfList.add( *ws.pdf(Form("%sTot_Jpsi_%s", pdfType.c_str(), "PP")) );  }
  if (incBkg)  { pdfList.add( *ws.pdf(Form("%sTot_Bkg_%s", pdfType.c_str(), "PP")) );   }
  if (!incJpsi && !incBkg) { cout << "[ERROR] User did not include any model, please fix your input settings!" << endl; return false; }
  RooAbsPdf *themodel = new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfList );
  ws.import(*themodel);
  ws.pdf(pdfName.c_str())->setNormRange("MassWindow");

  cout << "Model imported" << endl;
  setFixedVarsToContantVars(ws);
  cout << "Var set constatn" << endl;

  // save the initial values of the model we've just created
  RooArgSet* params = (RooArgSet*) themodel->getParameters(RooArgSet(*ws.var("invMass")));
  ws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);
  
  return true;
};

bool addBackgroundMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni) 
{
  cout << Form("[INFO] Implementing %s Background Mass Model", object.c_str()) << endl;

  // Import the Yield parameter
  if (!ws.var(Form("N_%s_%s", object.c_str(), "PP"))) { ws.factory(parIni[Form("N_%s_%s", object.c_str(), "PP")].c_str()); cout << "!!!!!!!!" << parIni[Form("N_%s_%s", object.c_str(), "PP")].c_str() << endl;}

  switch(model) 
    {  
    case (MassModel::Uniform): 

      // create the PDF           
      ws.factory(Form("Uniform::%s(%s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass"));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background Uniform PDF in %s included", object.c_str(), "PP") << endl; break;
 
    case (MassModel::Chebychev1): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background First Order Chebychev in %s", object.c_str(), "PP") << endl; return false;
      } 

      // create the variables for this model
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF           
      ws.factory(Form("Chebychev::%s(%s, {%s})", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 1st Order Chebychev PDF in %s included", object.c_str(), "PP") << endl; break;
 
    case (MassModel::Chebychev2): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Second Order Chebychev in %s", object.c_str(), "PP") << endl; return false;
      } 

      // create the variables for this model
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF           
      ws.factory(Form("Chebychev::%s(%s, {%s, %s})", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 2nd Order Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 

    case (MassModel::Chebychev3): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Third Order Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model 
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF                 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP"), 
		      Form("lambda3_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 3rd Order Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 

    case (MassModel::Chebychev4): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda4_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fourth Order Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF                 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP"), 
		      Form("lambda3_%s_%s", object.c_str(), "PP"), 
		      Form("lambda4_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 4th Order Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 

    case (MassModel::Chebychev5): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda4_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda5_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fifth Order Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF                 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP"), 
		      Form("lambda3_%s_%s", object.c_str(), "PP"), 
		      Form("lambda4_%s_%s", object.c_str(), "PP"), 
		      Form("lambda5_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 5th Order Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 

    case (MassModel::Chebychev6): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda4_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda5_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("lambda6_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Sixth Order Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda6_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF                 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP"), 
		      Form("lambda3_%s_%s", object.c_str(), "PP"), 
		      Form("lambda4_%s_%s", object.c_str(), "PP"), 
		      Form("lambda5_%s_%s", object.c_str(), "PP"), 
		      Form("lambda6_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 6th Order Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 

    case (MassModel::ExpChebychev1): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background First Order Exponential Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni["invMassNorm"].c_str() );
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF                 
      ws.factory(Form("RooFormulaVar::%s('@1*(@0) + 1.0', {%s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP"), "invMassNorm", 
		      Form("lambda1_%s_%s", object.c_str(), "PP")
		      ));               

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP")
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 1st Order Exponential Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 
   
    case (MassModel::ExpChebychev2): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP"))  
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Second Order Exponential Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model    
      ws.factory( parIni["invMassNorm"].c_str() );    
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF            
      ws.factory(Form("RooFormulaVar::%s('@2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP"), "invMassNorm",
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP")
		      ));             

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP")
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 2nd Order Exponential Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 
   
    case (MassModel::ExpChebychev3): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda3_%s_%s", object.c_str(), "PP"))  
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Third Order Exponential Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model  
      ws.factory( parIni["invMassNorm"].c_str() );      
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF
      ws.factory(Form("RooFormulaVar::%s('@3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP"), "invMassNorm",
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP"), 
		      Form("lambda3_%s_%s", object.c_str(), "PP")
		      ));        

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP")
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 3rd Order Exponential Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 
   
    case (MassModel::ExpChebychev4): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda3_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda4_%s_%s", object.c_str(), "PP"))    
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fourth Order Exponential Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model   
      ws.factory( parIni["invMassNorm"].c_str() );     
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF                 
      ws.factory(Form("RooFormulaVar::%s('@4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP"), "invMassNorm",
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP"), 
		      Form("lambda3_%s_%s", object.c_str(), "PP"), 
		      Form("lambda4_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP")
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 4th Order Exponential Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 
   
    case (MassModel::ExpChebychev5): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda3_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda4_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda5_%s_%s", object.c_str(), "PP"))    
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fifth Order Exponential Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model     
      ws.factory( parIni["invMassNorm"].c_str() );   
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF       
      ws.factory(Form("RooFormulaVar::%s('@5*(16.0*@0*@0*@0*@0*@0 - 20.0*@0*@0*@0 + 5.0*@0) + @4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP"), "invMassNorm", 
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP"), 
		      Form("lambda3_%s_%s", object.c_str(), "PP"), 
		      Form("lambda4_%s_%s", object.c_str(), "PP"), 
		      Form("lambda5_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP")
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 5th Order Exponential Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 
   
    case (MassModel::ExpChebychev6): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda3_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda4_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda5_%s_%s", object.c_str(), "PP")),
            parIni.count(Form("lambda6_%s_%s", object.c_str(), "PP"))     
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Sixth Order Exponential Chebychev in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model    
      ws.factory( parIni["invMassNorm"].c_str() );    
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("lambda6_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF
      ws.factory(Form("RooFormulaVar::%s('@6*(32.0*@0*@0*@0*@0*@0*@0 - 48.0*@0*@0*@0*@0 + 18.0*@0*@0 - 1.0) + @5*(16.0*@0*@0*@0*@0*@0 - 20.0*@0*@0*@0 + 5.0*@0) + @4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP"), "invMassNorm",
		      Form("lambda1_%s_%s", object.c_str(), "PP"), 
		      Form("lambda2_%s_%s", object.c_str(), "PP"), 
		      Form("lambda3_%s_%s", object.c_str(), "PP"), 
		      Form("lambda4_%s_%s", object.c_str(), "PP"), 
		      Form("lambda5_%s_%s", object.c_str(), "PP"), 
		      Form("lambda6_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASSPol_%s_%s", object.c_str(), "PP")
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background 6th Order Exponential Chebychev PDF in %s included", object.c_str(), "PP") << endl; break; 

    case (MassModel::Exponential): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Exponential in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF                 
      ws.factory(Form("Exponential::%s(%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Background Exponential PDF in %s included", object.c_str(), "PP") << endl; break;
   
    case (MassModel::VWG):  

      gROOT->ProcessLine(".L ./Macros/Utilities/VWGPdf.cxx+");
      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("a_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("b_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("c_%s_%s", object.c_str(), "PP"))
             )) {
	cout << Form("[ERROR] Initial parameters where not found for %s VWG Model in %s", object.c_str(), "PP") << endl; return false; 
      }
      // create the variables for this model
      ws.factory( parIni["invMassNorm"].c_str() );    
      ws.factory( parIni[Form("a_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("b_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("c_%s_%s", object.c_str(), "PP")].c_str() );

      ws.factory(Form("VWGPdf::%s(%s, %s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass",
                      Form("a_%s_%s", object.c_str(), "PP"),
                      Form("b_%s_%s", object.c_str(), "PP"),
                      Form("c_%s_%s", object.c_str(), "PP")
                      ));
	cout << Form("VWGPdf::%s(%s, %s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass",
                      Form("a_%s_%s", object.c_str(), "PP"),
                      Form("b_%s_%s", object.c_str(), "PP"),
                      Form("c_%s_%s", object.c_str(), "PP")
                      ) << endl;


      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

	cout << Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ) << endl;

      cout << Form("[INFO] %s VWG PDF in %s included", object.c_str(), "PP") << endl; break;
      
    default :
      
      cout << "[ERROR] Selected Background Mass Model: " << parIni[Form("Model_%s_%s", object.c_str(), "PP")] << " has not been implemented" << endl; return false;
    
    }
  
  return true;
};


bool addSignalMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni) 
{
  cout << Form("[INFO] Implementing %s Mass Model", object.c_str()) << endl;
  
  std::string lb = Form("_%s_%s", object.c_str(), "PP");
  RooAbsPdf* pdf = NULL;
  std::string s1="";
  std::string s2="";

  // Import the Yield parameter
  if (!ws.var(Form("N_%s_%s", object.c_str(), "PP"))) { ws.factory(parIni[Form("N_%s_%s", object.c_str(), "PP")].c_str()); cout << "!!!!!!!" << parIni[Form("N_%s_%s", object.c_str(), "PP")].c_str() << endl;}

  switch(model) 
    {    
    case (MassModel::SingleGaussian): 
      
      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Single Gaussian Model in %s", object.c_str(), "PP") << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF                       
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass", 
		      Form("m_%s_%s", object.c_str(), "PP"), 
		      Form("sigma1_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));
    
      cout << Form("[INFO] %s Single Gaussian PDF in %s included", object.c_str(), "PP") << endl; break;  
      
    case (MassModel::DoubleGaussian): 

      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("sigma2_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("f_%s_%s", object.c_str(), "PP")) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Double Gaussian Model in %s", object.c_str(), "PP") << endl; return false; 
      }

      // create the variables for this model              
      ws.factory( parIni[Form("f_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), "PP")].c_str() ); 
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("sigma2_%s_%s", object.c_str(), "PP")].c_str() );

      // create the two PDFs             
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS1_%s_%s", object.c_str(), "PP"), "invMass", 
                      Form("m_%s_%s", object.c_str(), "PP"), 
                      Form("sigma1_%s_%s", object.c_str(), "PP")
                      ));
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS2_%s_%s", object.c_str(), "PP"), "invMass", 
                      Form("m_%s_%s", object.c_str(), "PP"), 
                      Form("sigma2_%s_%s", object.c_str(), "PP")
                      ));

      // Sum the PDFs to get the signal PDF
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"),
		      Form("f_%s_%s", object.c_str(), "PP"),
		      Form("pdfMASS1_%s_%s", object.c_str(), "PP"),
		      Form("pdfMASS2_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Double Gaussian PDF in %s included", object.c_str(), "PP") << endl; break; 

    case (MassModel::SingleCrystalBall):  

      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("alpha_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("n_%s_%s", object.c_str(), "PP"))
             )) {
	cout << Form("[ERROR] Initial parameters where not found for %s Single Crystal Ball Model in %s", object.c_str(), "PP") << endl; return false; 
      }

      // create the variables for this model             
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), "PP")].c_str() );

      // create the PDF              
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass", 
		      Form("m_%s_%s", object.c_str(), "PP"), 
		      Form("sigma1_%s_%s", object.c_str(), "PP"),
		      Form("alpha_%s_%s", object.c_str(), "PP"),
		      Form("n_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Single Crystal Ball PDF in %s included", object.c_str(), "PP") << endl; break;
      
    case (MassModel::ExtendedCrystalBall):  

      gROOT->ProcessLine(".L ./Macros/Utilities/RooExtCBShape.cxx+");
      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("sigma1_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("alpha_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("n_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("alpha2_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("n2_%s_%s", object.c_str(), "PP"))
             )) {
	cout << Form("[ERROR] Initial parameters where not found for %s Extended Crystal Ball Model in %s", object.c_str(), "PP") << endl; return false; 
      }
	cout << "Defined var for ext : " << Form("alpha2_%s_%s", object.c_str(), "PP") << endl;
	cout << "NB in parIni : " << parIni.count(Form("alpha2_%s_%s", object.c_str(), "PP")) << endl;
	cout << "Value in ParIni : " << parIni[Form("alpha2_%s_%s", object.c_str(), "PP")].c_str() << endl;
	cout << "Value in ParIni : " << parIni[Form("m_%s_%s", object.c_str(), "PP")].c_str() << endl;
	cout << "Value in ParIni : " << parIni[Form("sigma1_%s_%s", object.c_str(), "PP")].c_str() << endl;
	cout << "Value in ParIni : " << parIni[Form("alpha_%s_%s", object.c_str(), "PP")].c_str() << endl;
	cout << "Value in ParIni : " << parIni[Form("n_%s_%s", object.c_str(), "PP")].c_str() << endl;
	cout << "Value in ParIni : " << parIni[Form("n2_%s_%s", object.c_str(), "PP")].c_str() << endl;

      // create the variables for this model
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), "PP")].c_str() );
      s1="alpha2_Jpsi_PP[2.0000,0.5000,30.0000]";
//      ws.factory(s1.c_str());
      ws.factory( parIni[Form("alpha2_%s_%s", object.c_str(), "PP")].c_str() );
      s2="n2_Jpsi_PP[1.8000,0.5000,10.0000]";
//      ws.factory(s2.c_str());
      ws.factory( parIni[Form("n2_%s_%s", object.c_str(), "PP")].c_str() );

cout << "Looking for var : " << ("alpha2"+lb).c_str() << endl;
cout << ("pdfMASS"+lb).c_str() << endl;

      // create the PDF
/*      pdf = new RooExtCBShape(("pdfMASS"+lb).c_str(), ("pdfMASS"+lb).c_str(),
                              *ws.var("invMass"),
                              *ws.var(("m"+lb).c_str()),
                              *ws.var(("sigma1"+lb).c_str()),
                              *ws.var(("alpha"+lb).c_str()),
                              *ws.var(("n"+lb).c_str()),
                              *ws.var(("alpha2"+lb).c_str()),
                              *ws.var(("n2"+lb).c_str())
                              );*/
cout << "pdf created" << endl;
//      if (pdf) { ws.import(*pdf); }
cout << "pdf imported" << endl;

      ws.factory(Form("RooExtCBShape::%s(%s, %s, %s, %s, %s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass",
                      Form("m_%s_%s", object.c_str(), "PP"),
                      Form("sigma1_%s_%s", object.c_str(), "PP"),
                      Form("alpha_%s_%s", object.c_str(), "PP"),
                      Form("n_%s_%s", object.c_str(), "PP"),
                      Form("alpha2_%s_%s", object.c_str(), "PP"),
                      Form("n2_%s_%s", object.c_str(), "PP")
                      ));

	cout << Form("RooExtCBShape::%s(%s, %s, %s, %s, %s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass",
                      Form("m_%s_%s", object.c_str(), "PP"),
                      Form("sigma1_%s_%s", object.c_str(), "PP"),
                      Form("alpha_%s_%s", object.c_str(), "PP"),
                      Form("n_%s_%s", object.c_str(), "PP"),
                      Form("alpha2_%s_%s", object.c_str(), "PP"),
                      Form("n2_%s_%s", object.c_str(), "PP")
                      ) << endl;

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

	cout << Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ) << endl;

      cout << Form("[INFO] %s Extended Crystal Ball PDF in %s included", object.c_str(), "PP") << endl; break;
      
    case (MassModel::DoubleCrystalBall): 
      
      // check that all input parameters are defined
        if (!(
              parIni.count(Form("m_%s_%s", object.c_str(), "PP")) &&
              parIni.count(Form("sigma1_%s_%s", object.c_str(), "PP")) &&
              parIni.count(Form("sigma2_%s_%s", object.c_str(), "PP")) &&
              parIni.count(Form("alpha_%s_%s", object.c_str(), "PP")) &&
              parIni.count(Form("alpha2_%s_%s", object.c_str(), "PP")) &&
              parIni.count(Form("n_%s_%s", object.c_str(), "PP")) &&
              parIni.count(Form("n2_%s_%s", object.c_str(), "PP")) &&
              parIni.count(Form("f_%s_%s", object.c_str(), "PP"))
              )) {
          cout << Form("[ERROR] Initial parameters where not found for %s Double Crystal Ball Model in %s", object.c_str(), "PP") << endl; return false;
        }
	cout << "Defined var for ext : " << Form("alpha2_%s_%s", object.c_str(), "PP") << endl;
	cout << "NB in parIni : " << parIni.count(Form("alpha2_%s_%s", object.c_str(), "PP")) << endl;
	cout << "Value in ParIni : " << parIni[Form("alpha2_%s_%s", object.c_str(), "PP")].c_str() << endl;
        
      // create the variables for this model             
      ws.factory( parIni[Form("f_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("alpha2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("n2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("sigma2_%s_%s", object.c_str(), "PP")].c_str() );

      // create the two PDFs
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS1_%s_%s", object.c_str(), "PP"), "invMass", 
                      Form("m_%s_%s", object.c_str(), "PP"), 
                      Form("sigma1_%s_%s", object.c_str(), "PP"),
                      Form("alpha_%s_%s", object.c_str(), "PP"),
                      Form("n_%s_%s", object.c_str(), "PP")
                      ));
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS2_%s_%s", object.c_str(), "PP"), "invMass", 
                      Form("m_%s_%s", object.c_str(), "PP"), 
                      Form("sigma2_%s_%s", object.c_str(), "PP"),
                      Form("alpha2_%s_%s", object.c_str(), "PP"),
                      Form("n2_%s_%s", object.c_str(), "PP")
                      ));

      // Sum the PDFs to get the signal PDF
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"),
		      Form("f_%s_%s", object.c_str(), "PP"),
		      Form("pdfMASS1_%s_%s", object.c_str(), "PP"),
		      Form("pdfMASS2_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Double Crystal Ball PDF in %s included", object.c_str(), "PP") << endl; break;
      
    case (MassModel::GaussianAndCrystalBall):

      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("sigma2_%s_%s", object.c_str(), "PP")) && 
            parIni.count(Form("alpha_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("n_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("f_%s_%s", object.c_str(), "PP"))  
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Gaussian and Crystal Ball Model in %s", object.c_str(), "PP") << endl; return false;
      } 

      // create the variables for this model             
      ws.factory( parIni[Form("f_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("sigma2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), "PP")].c_str() );
      
      // create the two PDFs
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS1_%s_%s", object.c_str(), "PP"), "invMass",
                        Form("m_%s_%s", object.c_str(), "PP"),
                        Form("sigma1_%s_%s", object.c_str(), "PP"),
                        Form("alpha_%s_%s", object.c_str(), "PP"),
                        Form("n_%s_%s", object.c_str(), "PP")
                        ));
        
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS2_%s_%s", object.c_str(), "PP"), "invMass",
                      Form("m_%s_%s", object.c_str(), "PP"), 
                      Form("sigma2_%s_%s", object.c_str(), "PP")
                      ));
 
      // Sum the PDFs to get the signal PDF 
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"),
		      Form("f_%s_%s", object.c_str(), "PP"),
		      Form("pdfMASS1_%s_%s", object.c_str(), "PP"),
		      Form("pdfMASS2_%s_%s", object.c_str(), "PP")
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s Gaussian and Crystal Ball PDF in %s included", object.c_str(), "PP") << endl; break;
	
    case (MassModel::NA60):  

      gROOT->ProcessLine(".L ./Macros/Utilities/NA60Shape.cxx+");
      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("sigma1_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("alpha_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("p1_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("p2_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("p3_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("alpha2_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("p12_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("p22_%s_%s", object.c_str(), "PP")) &&
            parIni.count(Form("p32_%s_%s", object.c_str(), "PP"))
             )) {
	cout << Form("[ERROR] Initial parameters where not found for %s NA60 Model in %s", object.c_str(), "PP") << endl; return false; 
      }
      // create the variables for this model
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("p1_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("p2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("p3_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("alpha2_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("p12_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("p22_%s_%s", object.c_str(), "PP")].c_str() );
      ws.factory( parIni[Form("p32_%s_%s", object.c_str(), "PP")].c_str() );

      ws.factory(Form("NA60Shape::%s(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), "PP"), "invMass",
                      Form("m_%s_%s", object.c_str(), "PP"),
                      Form("sigma1_%s_%s", object.c_str(), "PP"),
                      Form("alpha_%s_%s", object.c_str(), "PP"),
                      Form("p1_%s_%s", object.c_str(), "PP"),
                      Form("p2_%s_%s", object.c_str(), "PP"),
                      Form("p3_%s_%s", object.c_str(), "PP"),
                      Form("alpha2_%s_%s", object.c_str(), "PP"),
                      Form("p12_%s_%s", object.c_str(), "PP"),
                      Form("p22_%s_%s", object.c_str(), "PP"),
                      Form("p32_%s_%s", object.c_str(), "PP")
                      ));


      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), "PP"),
                      Form("pdfMASS_%s_%s", object.c_str(), "PP"),
                      Form("N_%s_%s", object.c_str(), "PP")
                      ));

      cout << Form("[INFO] %s NA60 PDF in %s included", object.c_str(), "PP") << endl; break;
      
    default :

      cout << "[ERROR] Selected Signal Mass Model: " << parIni[Form("Model_%s_%s", object.c_str(), "PP")] << " has not been implemented" << endl; return false;

    }
  
  return true;
};


void setMassDefaultParameters(map<string, string> &parIni, double numEntries)
{

  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;

  // DEFAULT SINGLE AND DOUBLE RATIO PARAMETERS

  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", "PP"))==0 || parIni[Form("N_Jpsi_%s", "PP")]=="") { 
    parIni[Form("N_Jpsi_%s", "PP")]  = Form("%s[%.10f,%.10f,%.10f]", Form("N_Jpsi_%s", "PP"), numEntries, -2.0*numEntries, 2.0*numEntries);
  }
  if (parIni.count(Form("N_Bkg_%s", "PP"))==0 || parIni[Form("N_Bkg_%s", "PP")]=="") { 
    parIni[Form("N_Bkg_%s", "PP")]  = Form("%s[%.10f,%.10f,%.10f]", Form("N_Bkg_%s", "PP"), numEntries, -2.0*numEntries, 2.0*numEntries);
  }

  // DEFAULT SIGNAL MASS MODEL PARAMETERS 
  if (parIni.count(Form("m_Jpsi_%s", "PP"))==0 || parIni[Form("m_Jpsi_%s", "PP")]=="") {
    parIni[Form("m_Jpsi_%s", "PP")] = Form("%s[%.6f,%.6f,%.6f]", Form("m_Jpsi_%s", "PP"), Mass.JPsi, Mass.JPsi-0.1, Mass.JPsi+0.1);
  }
  if (parIni.count(Form("sigma1_Jpsi_%s", "PP"))==0 || parIni[Form("sigma1_Jpsi_%s", "PP")]=="") {
    parIni[Form("sigma1_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma1_Jpsi_%s", "PP"), 0.05, 0.005, 0.2);
  }
  if (parIni.count(Form("rSigma21_Jpsi_%s", "PP"))==0) {
    if (parIni.count(Form("sigma2_Jpsi_%s", "PP"))==0 || parIni[Form("sigma2_Jpsi_%s", "PP")]=="") {
      parIni[Form("sigma2_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma2_Jpsi_%s", "PP"), 0.04, 0.01, 0.10);
    }
  } else {
    if (parIni[Form("rSigma21_Jpsi_%s", "PP")]=="") {
      parIni[Form("rSigma21_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("rSigma21_Jpsi_%s", "PP"), 2.0, 1.0, 4.0);
    }
    parIni[Form("sigma2_Jpsi_%s", "PP")] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("sigma2_Jpsi_%s", "PP"),
                                                                parIni[Form("rSigma21_Jpsi_%s", "PP")].c_str(), Form("sigma1_Jpsi_%s", "PP" ));
  }      
  if (parIni.count(Form("alpha_Jpsi_%s", "PP"))==0 || parIni[Form("alpha_Jpsi_%s", "PP")]=="") {
    parIni[Form("alpha_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha_Jpsi_%s", "PP"), 2.0, 0.5, 30.0);
//    parIni[Form("alpha_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha_Jpsi_%s", "PP"), -0.4, -1.0, 30.0);
  }
  if (parIni.count(Form("alpha2_Jpsi_%s", "PP"))==0) {
	cout << "alpha2 @0 case" << endl;
    parIni[Form("alpha2_Jpsi_%s", "PP")] = Form("RooFormulaVar::%s('@0',{%s})", Form("alpha2_Jpsi_%s", "PP"), Form("alpha_Jpsi_%s", "PP"));
  }
  else if (parIni[Form("alpha2_Jpsi_%s", "PP")]=="")
  {
	cout << "alpha2 other" << endl;
    parIni[Form("alpha2_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha2_Jpsi_%s", "PP"), 2.0, 0.5, 30.0);
//    parIni[Form("alpha2_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha2_Jpsi_%s", "PP"), 2.0, -1.0, 30.0);
  }else{
	cout << "alpha2 else" << endl;
  }
  if (parIni.count(Form("n_Jpsi_%s", "PP"))==0 || parIni[Form("n_Jpsi_%s", "PP")]=="") {
    parIni[Form("n_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("n_Jpsi_%s", "PP"), 1.8, 0.5, 10.0);
  }
  if (parIni.count(Form("n2_Jpsi_%s", "PP"))==0) {
    parIni[Form("n2_Jpsi_%s", "PP")] = Form("RooFormulaVar::%s('@0',{%s})", Form("n2_Jpsi_%s", "PP"), Form("n_Jpsi_%s", "PP"));
  }
  else if (parIni[Form("n2_Jpsi_%s", "PP")]=="")
  {
    parIni[Form("n2_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("n2_Jpsi_%s", "PP"), 12.0, 0.5, 100.0);
  }
  if (parIni.count(Form("f_Jpsi_%s", "PP"))==0 || parIni[Form("f_Jpsi_%s", "PP")]=="") {
    parIni[Form("f_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("f_Jpsi_%s", "PP"), 0.5, 0.0, 1.0);
  }


  if (parIni.count(Form("p1_Jpsi_%s", "PP"))==0 || parIni[Form("p1_Jpsi_%s", "PP")]=="") {
    parIni[Form("p1_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("p1_Jpsi_%s", "PP"), 0.23, 0.01, 10.0);
  }
  if (parIni.count(Form("p2_Jpsi_%s", "PP"))==0 || parIni[Form("p2_Jpsi_%s", "PP")]=="") {
    parIni[Form("p2_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("p2_Jpsi_%s", "PP"), 1.2, 0.01, 10.0);
  }
  if (parIni.count(Form("p3_Jpsi_%s", "PP"))==0 || parIni[Form("p3_Jpsi_%s", "PP")]=="") {
    parIni[Form("p3_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("p3_Jpsi_%s", "PP"), 0.04, 0.01, 10.0);
  }

  if (parIni.count(Form("a_Bkg_%s", "PP"))==0 || parIni[Form("a_Bkg_%s", "PP")]=="") {
    parIni[Form("a_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("a_Bkg_%s", "PP"), 1., -10., 10.0);
  }
  if (parIni.count(Form("b_Bkg_%s", "PP"))==0 || parIni[Form("b_Bkg_%s", "PP")]=="") {
    parIni[Form("b_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("b_Bkg_%s", "PP"), 1., 0.0, 100.0);
  }
  if (parIni.count(Form("c_Bkg_%s", "PP"))==0 || parIni[Form("c_Bkg_%s", "PP")]=="") {
    parIni[Form("c_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("c_Bkg_%s", "PP"), 1., -100., 100.0);
  }

  if (parIni.count(Form("p12_Jpsi_%s", "PP"))==0 || parIni[Form("p12_Jpsi_%s", "PP")]=="") {
    parIni[Form("p12_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("p12_Jpsi_%s", "PP"), 0.18, 0.01, 10.0);
  }
  if (parIni.count(Form("p22_Jpsi_%s", "PP"))==0 || parIni[Form("p22_Jpsi_%s", "PP")]=="") {
    parIni[Form("p22_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("p22_Jpsi_%s", "PP"), 1.3, 0.01, 10.0);
  }
  if (parIni.count(Form("p32_Jpsi_%s", "PP"))==0 || parIni[Form("p32_Jpsi_%s", "PP")]=="") {
    parIni[Form("p32_Jpsi_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("p32_Jpsi_%s", "PP"), 0.06, 0.01, 10.0);
  }

  // DEFAULT BACKGROUND MASS MODEL PARAMETERS
  if (parIni[Form("Model_Bkg_%s","PP")].find("ExpChebychev")!=std::string::npos) {
    if (parIni.count(Form("lambda1_Bkg_%s", "PP"))==0 || parIni[Form("lambda1_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda1_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", "PP"), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda2_Bkg_%s", "PP"))==0 || parIni[Form("lambda2_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda2_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda2_Bkg_%s", "PP"), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda3_Bkg_%s", "PP"))==0 || parIni[Form("lambda3_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda3_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda3_Bkg_%s", "PP"), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda4_Bkg_%s", "PP"))==0 || parIni[Form("lambda4_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda4_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda4_Bkg_%s", "PP"), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda5_Bkg_%s", "PP"))==0 || parIni[Form("lambda5_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda5_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda5_Bkg_%s", "PP"), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda6_Bkg_%s", "PP"))==0 || parIni[Form("lambda6_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda6_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda6_Bkg_%s", "PP"), 0.0, -10.0, 10.0);
    }
  }
  else if (parIni[Form("Model_Bkg_%s","PP")].find("Chebychev")!=std::string::npos) {
    if (parIni.count(Form("lambda1_Bkg_%s", "PP"))==0 || parIni[Form("lambda1_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda1_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", "PP"), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda2_Bkg_%s", "PP"))==0 || parIni[Form("lambda2_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda2_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda2_Bkg_%s", "PP"), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda3_Bkg_%s", "PP"))==0 || parIni[Form("lambda3_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda3_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda3_Bkg_%s", "PP"), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda4_Bkg_%s", "PP"))==0 || parIni[Form("lambda4_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda4_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda4_Bkg_%s", "PP"), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda5_Bkg_%s", "PP"))==0 || parIni[Form("lambda5_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda5_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda5_Bkg_%s", "PP"), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda6_Bkg_%s", "PP"))==0 || parIni[Form("lambda6_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda6_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda6_Bkg_%s", "PP"), 0.0, -2.0, 2.0);
    }
  } 
  else if (parIni[Form("Model_Bkg_%s","PP")].find("Exponential")!=std::string::npos) {
    if (parIni.count(Form("lambda1_Bkg_%s", "PP"))==0 || parIni[Form("lambda1_Bkg_%s", "PP")]=="") { 
      parIni[Form("lambda1_Bkg_%s", "PP")] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", "PP"), 0.0, -100.0, 100.0);
    }
  }
 
};


#endif // #ifndef buildCharmoniaMassModel_C
