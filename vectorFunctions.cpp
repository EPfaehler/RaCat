//template<typename T>
void getGreaterElements(vector<float> &vectorOfMatrixElem, float &value){
    typename vector<float>::iterator greaterThan;
    greaterThan = remove_if(vectorOfMatrixElem.begin(), vectorOfMatrixElem.end(), bind2nd(greater<float>(), value));
    vectorOfMatrixElem.erase(greaterThan, vectorOfMatrixElem.end());
}
