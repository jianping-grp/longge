from . import models
from dynamic_rest import serializers
from rest_framework.permissions import AllowAny

class SubmitParamsSerializer(serializers.DynamicModelSerializer):
    class Meta:
        model = models.SubmitParamter
        exclude = []
