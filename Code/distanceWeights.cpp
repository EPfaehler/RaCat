
float calculateManhattanNorm2D(int directionX, int directionY, vector<double> spacing) {
    float actDistance = abs(spacing[0]*directionX) + abs(spacing[1]*directionY);
    return actDistance;
}


float calculateEuclidianNorm2D(int directionX, int directionY, vector<double> spacing) {
    float actDistance = sqrt(pow(directionX*spacing[0],2) + pow(directionY*spacing[1],2));
    return actDistance;
}

float calculateManhattanNorm3D(int directionX, int directionY, int directionZ, vector<double> spacing) {
	float actDistance = abs(spacing[0] * directionX) + abs(spacing[1] * directionY) + abs(spacing[2] * directionZ);
	return actDistance;
}


float calculateEuclidianNorm3D(int directionX, int directionY, int directionZ, vector<double> spacing) {
	float actDistance = sqrt(pow(directionX*spacing[0], 2) + pow(directionY*spacing[1], 2) + pow(directionZ*spacing[2], 2));
	return actDistance;
}
