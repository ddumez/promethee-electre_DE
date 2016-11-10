/* ---------------------- -*- Decision Engineering -*- ----------------------- */
//                                 main.cpp                                    //
//                                   C++                                       //
//                                  Projet                                     //
/* --------------------------------------------------------------------------- */
//              Réalisation de deux méthodes ELECTRE et PROMETHEE              //
/* --------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>

using namespace std;

#define NBVAR 200 //entre 1 et 200
#define NBCRIT 10 //entre 1 et 10

//------------------ fonctions generiques ---------------------
/**
* \brief parse le fichier
*
* \param[in] nom le nom du fichier de donnee
* \param[out] data le tableau dans lequel on met les donnees
*/
void lecture (const char* const nom, double data[NBVAR][NBCRIT]);
/**
* \brief affiche les valeures selon les differents critere de toutes les possibilites
*
* \param[in] data le tableau qui contien les valeurs des fonction eco des possibilite
*/
void affichage (const double data[NBVAR][NBCRIT]);

//------------------ fonction electre -------------------------
/**
* \brief calcule les indices de concordance et de discordance
*
* \param[in] data les valeurs des fonctions objectif des differentes possibilite
* \param[out] c le tableau des indices de concordance
* \param[out] d le tableau des indices de discordance
*/
void calcul_discordance_concordance(const double data[NBVAR][NBCRIT], double c[NBVAR][NBVAR], double d[NBVAR][NBVAR]);
/**
* \brief determine pour chaque couple si aj surclasse ak
*
* \param[in] data les donnes sur les possibilitees
* \param[out] S la matrice de surclassement
* \param[in] seuilc seuil a partir duqel la concordance est significative
* \param[in] seuild seuil a partir duqel la discordance est significative
*/
void caclul_electre(const double data[NBVAR][NBCRIT], bool S[NBVAR][NBVAR], const double seuilc, const double seuild);
/**
* \brief trie le tableau possibilite en fonction du nombre de fois que la possibilite est non surclasse
*
* \param[in] non_surclasse nombre de fois que chaque variable est non surclasse
* \param[out] possibilite indice des possibilite dans le tableau non_surclasse
*/
void trie_non_surclasse(const int non_surclasse[NBVAR], int possibilite[NBVAR]);

//------------------ fonctions promethee ----------------------
/**
* \brief caclcule la valeur de la fonction de preference normlise de promethee
*
* \param[in] diff la difference entre les deux valeures a comparer
* \param[in] p le seuil minimal de significativite
* \param[in] q la borne max de significativite
* \return la valeur de la fonction de prefence de promethee
*/
float fonction_de_preference (const double diff, const int p, const int q);
/**
* \brief calcule la preference pi pour toutes les paires de possibilite
*
* \param[in] data les valeurs de toutes les preferences
* \param[out] preference la preference pi pour toutes les paires de possibilite
*/
void calul_preference(const double data[NBVAR][NBCRIT], double prefence[NBVAR][NBVAR], const int p, const int q);
/**
* \brief calcul les flots par rapport aux preferences
*
* \param[in] preference les preference entre paire de possibilite
* \param[out] flotsortant phi majuscule +
* \param[out] flotentrant phi majuscule -
*/
void calcul_flots(const double prefence[NBVAR][NBVAR], double flotsortant[NBVAR], double flotentrant[NBVAR]);
/**
* \brief trie les possiblite en fonction de leur flot
*
* \param[in] flotsortant flots sortant des possibilites
* \param[in] flotentrant flots entrant des possibilites
* \param[out] possibilite la permutation des possiblites
*/
void trie_flot_net(const double flotsortant[NBVAR], const double flotentrant[NBVAR], int possibilite[NBVAR]);



int main(const int argc, const char* const argv[]) {
    int i,j,k,l; //compteur de boucle

    // Matrice de taille 100*10 destinée à stocker la matrice du fichier d'entrée
    double data[NBVAR][NBCRIT];
    
    // Indices de concordance et discordance pour ELECTRE
    const double concordance[4] = {0.6,0.7,0.8,0.9};
    const double discordance[5] = {0.1,0.15,0.2,0.25,0.3};
//    const double concordance[1] = {0.7};
//    const double discordance[1] = {0.15};
 

    // Seuils d'indifférence pour PROMETHEE
    const int p[3] = {10,15,20};
    const int q[3] = {90,85,80};
    
    // Lecture des données    
    lecture(argv[1], data);

    // Vérification de la lecture
    affichage(data);
    
    // ELECTRE
	    cout<<"\n\ncalul ELECTRE : "<<endl;
	    int non_surclasse[NBVAR]; //compte le nombre de fois que la variable est non surclasse
	    for(i = 0; i<NBVAR; ++i) {non_surclasse[i] = 0;}
		bool S[NBVAR][NBVAR]; //indique si aiSak
		int possibilite[NBVAR]; //pour ensuite calsser les variables par rapoort au nombre de fois qu'elles sont non-surclasse
		for(i = 0; i<NBVAR; ++i) {possibilite[i] = i;} //initialise a la permutation identite

		for(i = 0; i<4; ++i) { //pour toutes le concordances
			for(j = 0; j<5; ++j) { //pour toutes les discordances
				//calcul des surclassement avec ces valeurs des bornes
				caclul_electre(data, S, concordance[i], discordance[j]);

				//recherche des variable non surclasse
				for(k = 0; k<NBVAR; ++k) {
					l = 0;
					while( (l<NBVAR) && ( (!S[l][k]) || (l == k) ) ) {
						//tant que l ne surclasse pas k
						++l;
					}
					if (NBVAR <= l) {++ non_surclasse[k];} //si on est alle jusqu'au bout c'est que rien ne surclasse k
				}
			}
		}

		//trie_non_surclasse(non_surclasse, possibilite);
		cout<<"nombre de fois que chaque variable est non surclasse : "<<endl;
		for(i = 0; i<NBVAR; ++i) {cout<<possibilite[i]<<" ( "<<non_surclasse[possibilite[i]]<<" )"<<endl;}

	//PROMETHEE
		cout<<"\n\n calcul PROMETHEE"<<endl;
		double prefence[NBVAR][NBVAR];
		double flotsortant[NBVAR];
		double flotentrant[NBVAR];
		for(i = 0; i<3; ++i) { //pour toutes les paires de valeurs de seuils d'indiference
			calul_preference(data, prefence, p[i], q[i]);
			calcul_flots(prefence, flotsortant, flotentrant);

			trie_flot_net(flotsortant, flotentrant, possibilite);
			cout<<"flots net pour les valeurs "<<p[i]<<" et "<<q[i]<<" : "<<endl;
			for(j = 0; j<NBVAR; ++j) {cout<<possibilite[j]<<" : "<<flotsortant[possibilite[j]] - flotentrant[possibilite[j]]<<"\n";}
			cout<<"\n"<<endl;
		}

    //fin
	return 0;
}

// Lecture des données à partir du fichier "nom"
void lecture (const char* const nom, double data[NBVAR][NBCRIT]) {
		
	ifstream f(nom);
	
	// Lecture du fichier
	for (int i = 0; i < NBVAR; i++)
		for (int j = 0; j < NBCRIT; j++)
			f >> data[i][j];
    f.close(); 
}

// Affichage des données sur console 
void affichage (const double data[NBVAR][NBCRIT]) {
	
	cout<<"affichage des données"<<endl;
	for (int i = 0; i < NBVAR; i++) {
		for (int j = 0; j < NBCRIT; j++)
			cout << data[i][j] << "\t";
		cout << endl;
	}
}

//calul de la discordance et concordance
void calcul_discordance_concordance(const double data[NBVAR][NBCRIT], double c[NBVAR][NBVAR], double d[NBVAR][NBVAR]) {
	//variables
		int i, j, k;
		bool d0; //indique si la valeur de la disordance est 0
		double max;

	//debut
		//calul des concordances et discordance de chaque paire
		for (i = 0; i<NBVAR; ++i) {
			for(k = 0; k<NBVAR; ++k) {
				c[i][k] = 0;
				d0 = true;
				max = data[k][0] - data[i][0];
				for(j = 0; j<NBCRIT; ++j) {
					//calul indice de concordance
					if (data[i][j] >= data[k][j]) {
						c[i][k] += (double)(1.0/(double)(NBCRIT)); //poid de chaque critere egal
					}

					//calcul preparatoire de l'indice de discordance
					d0 = d0 && (data[i][j] >= data[k][j]); //si i domine k
					if (data[k][j] - data[i][j] > max) {max = data[k][j] - data[i][j];} //recherche de la plus grande difference
				}
				//*calcul indice de discordance
				if (d0) {
					d[i][k] = 0.0;
				} else {
					d[i][k] = (double)(max/98.0);
				}
			}
		}
	//fin
}
//calcul du surclassement selon electre
void caclul_electre(const double data[NBVAR][NBCRIT], bool S[NBVAR][NBVAR], const double seuilc, const double seuild) {
	//variable
		double c[NBVAR][NBVAR]; //tableau des indices de concordance
    	double d[NBVAR][NBVAR]; //tableau des indices de discordance
    	int i,k;
    
	//debut
    	//calcul des tableau de concordance et de discordance
    	calcul_discordance_concordance(data, c, d);

    	//interpretation
		for(i = 0; i < NBVAR; ++i) {
			for(k = 0; k < NBVAR; ++k) {
				S[i][k] = ( (c[i][k] >= seuilc) && (d[i][k] <= seuild) );
			}
		}
	//fin
}
//trie des possibilite en fonction du nombre de fois qu'elles sont non surclasse
void trie_non_surclasse(const int non_surclasse[NBVAR], int possibilite[NBVAR]) {
	int tmp;
	for(int i = 0; i<NBVAR; ++i) {
		for(int j = 0; j<NBVAR-1; ++j) {
			if (non_surclasse[possibilite[j]] < non_surclasse[possibilite[j+1]] ) {
				tmp = possibilite[j];
				possibilite[j] = possibilite[j+1];
				possibilite[j+1] = tmp;
			}
		}
	}
}


// Fonction de préférence
float fonction_de_preference (const double diff, const int p, const int q) {
	if (diff <= p) 
		return 0.0;
	else if (diff >= q)
		return 1.0;
	else 
	 	return 1.0*(diff - p)/(q - p);
}
//calcul de la preference globale
void calul_preference(const double data[NBVAR][NBCRIT], double prefence[NBVAR][NBVAR], const int p, const int q) {
	//variable
		int i,j,k;
	//debut
		for (i = 0; i<NBVAR; ++i) { //pour toutes les paires de possibilites
			for(j = 0; j<NBVAR; ++j) {
				prefence[i][j] = 0;
				for(k = 0; k<NBCRIT; ++k) { //on somme les pi_k
					prefence[i][j] += (1.0/NBCRIT) * fonction_de_preference(data[i][k] - data[j][k], p, q);
				}
			}
		}
	//fin
}
//calcul des flots
void calcul_flots(const double prefence[NBVAR][NBVAR], double flotsortant[NBVAR], double flotentrant[NBVAR]) {
	//variable
		int i,j;
	//debut
		for(i = 0; i<NBVAR; ++i) {
			flotsortant[i] = 0;
			flotentrant[i] = 0;
			for(j = 0; j<NBVAR; ++j) {
				flotsortant[i] += (1.0/(NBVAR-1)) * prefence[i][j];
				flotentrant[i] += (1.0/(NBVAR-1)) * prefence[j][i];
			}
		}
	//fin
}
//tri des possibilites
void trie_flot_net(const double flotsortant[NBVAR], const double flotentrant[NBVAR], int possibilite[NBVAR]) {
	int tmp;
	for(int i = 0; i<NBVAR; ++i) {
		for(int j = 0; j<NBVAR-1; ++j) {
			if (flotsortant[possibilite[j]] - flotentrant[possibilite[j]] < flotsortant[possibilite[j+1]] - flotentrant[possibilite[j+1]] ) {
				tmp = possibilite[j];
				possibilite[j] = possibilite[j+1];
				possibilite[j+1] = tmp;
			}
		}
	}
}